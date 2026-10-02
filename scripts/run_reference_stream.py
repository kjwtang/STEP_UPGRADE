#!/usr/bin/env python3
"""Restartable, bounded rain-only execution of an explicitly named preset.

No tracking or identification rules are changed. Only filename-hourly RAINNC
and projection-attribute grids are supported here. Graphs and checkpoints are
published together per chunk; labels are optional diagnostics, not the default.
"""
import argparse
from collections import Counter
from contextlib import contextmanager
import fcntl
import hashlib
from importlib.metadata import version, PackageNotFoundError
import json
import os
from pathlib import Path
import time
from types import SimpleNamespace
import uuid

import numpy as np

from step.identification import identify
from step.parallel_identification import FrameIdentifier
from step.tracking import load_tracking_state, save_tracking_state, track_with_graph
from frame_progress import FrameProgress
from run_npy_validation import disk, write_graph
from wrf_rain_input import manifest, read_wrf

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PRESET = ROOT / 'configs/rainfall_lineage_4km_1h_reference_v1.json'
SOURCE_FILES = ('step/tracking.py', 'step/identification.py', 'step/parallel_identification.py',
                'scripts/wrf_rain_input.py', 'scripts/run_reference_stream.py',
                'scripts/run_npy_validation.py', 'scripts/frame_progress.py')


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def atomic_json(path, value):
    temporary = path.with_name(f'.{path.name}.{uuid.uuid4().hex}.tmp')
    with temporary.open('w') as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def fingerprint(entries):
    records = []
    for stamp, path in entries:
        stat = path.stat()
        records.append(dict(timestamp=str(stamp), path=str(path.resolve()),
                            size=stat.st_size, mtime_ns=stat.st_mtime_ns))
    return records


def grid_contract(meta):
    return {key: meta[key] for key in
            ('source_shape', 'crop_origin_yx', 'dt_hours', 'normalized_units',
             'grid_attributes', 'time_source', 'grid_source')}


def payload_hashes(part):
    return {path.name: digest(path) for path in sorted(part.iterdir())
            if path.is_file() and path.name != 'COMMITTED.json'}


def committed_parts(out, contract):
    """Validate published chunks, including their output bytes, before resume."""
    parts = []
    offset = 0
    previous_last_time = -1
    expected_contract = contract_hash(contract)
    for path in sorted(out.glob('chunk_*')):
        record = json.loads((path / 'COMMITTED.json').read_text())
        if path.name != f'chunk_{offset:06d}' or record['offset'] != offset:
            raise ValueError('Published chunks are not contiguous')
        if record['frames'] < 1 or offset + record['frames'] > contract['hours']:
            raise ValueError('Invalid published chunk length')
        if record['contract_sha256'] != expected_contract:
            raise ValueError('Chunk execution contract differs')
        if record['files'] != payload_hashes(path):
            raise ValueError(f'Published output changed: {path.name}')
        if record['state_last_time'] != offset + record['frames'] - 1:
            raise ValueError('Checkpoint time does not match committed frames')
        if record['state_sequence_id'] != contract['sequence_id'] or record['state_last_time'] <= previous_last_time:
            raise ValueError('Checkpoint sequence or chronological order differs')
        previous_last_time = record['state_last_time']
        offset += record['frames']
        parts.append((path, record))
    return parts, offset


def contract_hash(contract):
    return hashlib.sha256(json.dumps(contract, sort_keys=True).encode()).hexdigest()


@contextmanager
def execution_lock(out):
    """Single-writer POSIX lock, automatically released after a killed process."""
    with (out / '.execution.lock').open('a') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError('Another process is executing this output directory') from error
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def dependency_versions():
    result = {}
    for name in ('numpy', 'scipy', 'xarray', 'netCDF4'):
        try:
            result[name] = version(name)
        except PackageNotFoundError:
            result[name] = 'not-installed'
    return result


def run(a):
    if a.hours < 1 or a.start < 0 or a.workers < 1 or a.chunk_frames < 1:
        raise ValueError('hours/chunk/workers must be positive; start must be nonnegative')
    if a.stop_after_chunks is not None and a.stop_after_chunks < 1:
        raise ValueError('stop-after-chunks must be positive')
    preset = json.loads(a.preset.read_text())
    if preset['input_contract'] != dict(grid_km=4, dt_hours=1,
                                       rain_units_after_normalization='mm/h'):
        raise ValueError('This runner requires the explicit 4-km, 1-hour, mm/h contract')
    entries = manifest(a.input, filename_time=True)
    selected = entries[a.start:a.start + a.hours + 1]
    if len(selected) != a.hours + 1:
        raise ValueError('N rain intervals need N+1 available cumulative snapshots')
    if any(b[0] - c[0] != np.timedelta64(1, 'h') for c, b in zip(selected, selected[1:])):
        raise ValueError('Planned snapshots must be consecutive hourly observations')
    reader = SimpleNamespace(input=a.input, wrf_time_source='filename',
        wrf_grid_source='attributes', inspect=False, variable=None,
        rain_kind='cumulative', start=a.start, hours=1, crop=0, units=None,
        dt_hours=1, negative_tolerance_mm=0)
    probe, metadata = read_wrf(reader, entries=entries)
    del probe
    grid = grid_contract(metadata)
    if any(grid['grid_attributes'][key] != 4000 for key in ('DX', 'DY')):
        raise ValueError('DX and DY must both be 4000 metres for this preset')
    # Total directory size is not a spatial contract; allow unrelated files to
    # arrive while requiring every planned snapshot to remain unchanged.
    grid['source_shape'] = grid['source_shape'][1:]
    contract = dict(format_version=1, preset=preset, sequence_id=a.sequence_id,
        start=a.start, hours=a.hours, save_labels=a.save_labels, grid=grid,
        snapshots=fingerprint(selected), dependency_versions=dependency_versions(),
        source_sha256={name: digest(ROOT / name) for name in SOURCE_FILES})
    out = a.output_dir
    if not a.resume:
        out.mkdir(parents=True, exist_ok=False)
    with execution_lock(out):
        return execute(a, contract, entries, reader, preset, grid)


def execute(a, contract, entries, reader, preset, grid):
    out = a.output_dir
    selected = entries[a.start:a.start + a.hours + 1]
    contract_digest = contract_hash(contract)
    if a.resume:
        saved = json.loads((out / 'execution_contract.json').read_text())
        if saved != contract:
            raise ValueError('Resume contract differs: preset, code, input, sequence or output mode changed')
        parts, offset = committed_parts(out, saved)
    else:
        atomic_json(out / 'execution_contract.json', contract)
        parts, offset = [], 0
    if (out / 'SUCCESS').exists():
        if offset != a.hours:
            raise ValueError('SUCCESS exists without all requested frames')
        marker = json.loads((out / 'SUCCESS').read_text())
        summary = json.loads((out / 'summary.json').read_text())
        if marker != dict(contract_sha256=contract_digest, frames=offset) or any((
                summary.get('contract_sha256') != contract_digest,
                summary.get('frames') != offset, summary.get('status') != 'completed',
                summary.get('sequence_id') != a.sequence_id,
                summary.get('nodes') != sum(r['nodes'] for _, r in parts),
                summary.get('edges') != sum(r['edges'] for _, r in parts))):
            raise ValueError('SUCCESS/summary does not match committed chunks')
        return summary
    state = load_tracking_state(parts[-1][0] / 'state.json') if parts else None
    if state is not None and (state.last_time != offset - 1 or state.sequence_id != a.sequence_id):
        raise ValueError('Latest checkpoint time or sequence differs')
    progress = FrameProgress(a.hours, enabled=not a.no_progress, initial_done=offset)
    progress.unit = 'hours'
    began = time.perf_counter()
    count = 0
    kwargs = dict(preset['tracking'])
    kwargs['km'] = kwargs.pop('max_displacement')
    identification = preset['identification']
    completed = False
    paused = False
    identifier = None
    try:
        if getattr(a, 'id_executor', 'fresh') == 'persistent':
            identifier = FrameIdentifier(grid['source_shape'], disk(identification['bridge_radius_cells']),
                capacity=a.chunk_frames, workers=a.workers,
                threshold=identification['threshold_mm_h'], min_size=identification['min_size_cells'])
        while offset < a.hours:
            reader.start = a.start + offset
            reader.hours = min(a.chunk_frames, a.hours - offset)
            size = reader.hours
            progress.update(offset, phase=f'reading frames {offset + 1}-{offset + size}')
            t = time.perf_counter()
            if fingerprint(entries[reader.start:reader.start + size + 1]) != contract['snapshots'][offset:offset + size + 1]:
                raise ValueError('Planned input changed before read')
            rain, meta = read_wrf(reader, entries=entries)
            current_grid = grid_contract(meta)
            current_grid['source_shape'] = current_grid['source_shape'][1:]
            if current_grid != grid:
                raise ValueError('Input grid/units/cadence contract changed between chunks')
            timings = {'read_seconds': time.perf_counter() - t}
            progress.update(phase='identification')
            t = time.perf_counter()
            if identifier is None:
                labels = identify(rain, disk(identification['bridge_radius_cells']),
                    threshold=identification['threshold_mm_h'],
                    min_size=identification['min_size_cells'], workers=a.workers,
                    parallel_mode='frames')
            else:
                labels = identifier.identify(rain)
            timings['identification_seconds'] = time.perf_counter() - t
            conflicts, decisions = [], []
            t = time.perf_counter()
            tracked, graph, state = track_with_graph(labels, rain, state=state,
                return_state=True, sequence_id=a.sequence_id,
                gap_conflict_callback=conflicts.extend, overlap_callback=decisions.extend,
                progress_callback=lambda done, total: progress.update(offset + done), **kwargs)
            timings['tracking_seconds'] = time.perf_counter() - t
            progress.update(phase='publishing chunk/checkpoint')
            t = time.perf_counter()
            # A killed writer leaves only a hidden staging directory. Resume
            # ignores it and never deletes user files or unpublished evidence.
            part = out / f'chunk_{offset:06d}'
            staging = out / f'.chunk_{offset:06d}.{uuid.uuid4().hex}.pending'
            staging.mkdir()
            write_graph(graph, staging)
            atomic_json(staging / 'metadata.json', meta)
            atomic_json(staging / 'gap_conflicts.json', conflicts)
            atomic_json(staging / 'adjacent_overlap.json', decisions)
            save_tracking_state(state, staging / 'state.json')
            if a.save_labels:
                np.save(staging / 'identified_labels.npy', labels)
                np.save(staging / 'tracked_labels.npy', tracked)
            record = dict(offset=offset, frames=size, nodes=len(graph.objects),
                workers=a.workers, id_executor=getattr(a, 'id_executor', 'fresh'),
                edges=len(graph.edges), gap_edges=sum(e.gap > 0 for e in graph.edges),
                events=dict(Counter(e.event for e in graph.events)),
                identified_sha256=hashlib.sha256(memoryview(labels)).hexdigest(),
                tracked_sha256=hashlib.sha256(memoryview(tracked)).hexdigest(),
                identified_frame_sha256=[hashlib.sha256(memoryview(f)).hexdigest() for f in labels],
                tracked_frame_sha256=[hashlib.sha256(memoryview(f)).hexdigest() for f in tracked],
                checkpoint_bytes=(staging / 'state.json').stat().st_size,
                state_last_time=state.last_time, state_sequence_id=state.sequence_id,
                historical_branches=len(state.family_parent),
                active=len(state.active), dormant=len(state.dormant),
                timings=timings, contract_sha256=contract_digest)
            record['files'] = payload_hashes(staging)
            timings['output_seconds'] = time.perf_counter() - t
            atomic_json(staging / 'COMMITTED.json', record)
            if part.exists():
                raise FileExistsError(part)
            os.replace(staging, part)
            parts.append((part, record))
            offset += size
            count += 1
            print(json.dumps({key: record[key] for key in
                ('offset', 'frames', 'nodes', 'edges', 'checkpoint_bytes', 'timings')}), flush=True)
            # Exercise actual serialization at every boundary, not merely an
            # in-memory object forwarded to the next call.
            state = load_tracking_state(part / 'state.json')
            del rain, labels, tracked, graph
            if a.stop_after_chunks is not None and count >= a.stop_after_chunks and offset < a.hours:
                atomic_json(out / 'PAUSED.json', dict(committed_frames=offset,
                    reason='requested stop-after-chunks; resume with unchanged execution contract'))
                progress.update(offset, phase='PAUSED after committed chunk; resume required')
                paused = True
                return dict(status='paused', committed_frames=offset)
        if fingerprint(selected) != contract['snapshots']:
            raise ValueError('Planned input changed during execution')
        totals = Counter()
        for _, record in parts:
            totals.update(record['timings'])
        summary = dict(status='completed', frames=offset, chunks=len(parts),
            nodes=sum(r['nodes'] for _, r in parts), edges=sum(r['edges'] for _, r in parts),
            gap_edges=sum(r['gap_edges'] for _, r in parts),
            events=dict(sum((Counter(r['events']) for _, r in parts), Counter())),
            stage_seconds=dict(totals), current_invocation_wall_seconds=time.perf_counter() - began,
            preset_id=preset['preset_id'], sequence_id=a.sequence_id,
            contract_sha256=contract_digest, last_time=state.last_time,
            historical_branches=len(state.family_parent),
            final_checkpoint_bytes=parts[-1][1]['checkpoint_bytes'],
            note='Filename cadence and attribute grid checked, not geographical or meteorological certification. Timings sum committed chunks across invocations; no whole-run equivalence inferred.')
        atomic_json(out / 'summary.json', summary)
        if (out / 'PAUSED.json').exists():
            atomic_json(out / 'PAUSED.json', dict(status='resolved', committed_frames=offset,
                                                reason='previous pause completed by resume'))
        atomic_json(out / 'SUCCESS', dict(contract_sha256=contract_digest, frames=offset))
        completed = True
        progress.finish(True)
        return summary
    finally:
        if identifier is not None:
            identifier.close(terminate=not (completed or paused))
        if not completed and not paused:
            progress.finish(False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input', type=Path)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--preset', type=Path, default=DEFAULT_PRESET)
    p.add_argument('--start', type=int, default=0, help='Preceding cumulative snapshot index')
    p.add_argument('--hours', type=int, required=True)
    p.add_argument('--sequence-id', required=True)
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--chunk-frames', type=int, default=24)
    p.add_argument('--id-executor', choices=['fresh', 'persistent'], default='fresh',
                   help='Opt-in reusable shared-mask pool; identical identification rules')
    p.add_argument('--save-labels', action='store_true')
    p.add_argument('--resume', action='store_true')
    p.add_argument('--stop-after-chunks', type=int)
    p.add_argument('--no-progress', action='store_true')
    print(json.dumps(run(p.parse_args()), indent=2), flush=True)


if __name__ == '__main__':
    main()
