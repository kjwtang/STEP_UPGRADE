"""Experimental tracking geometries; keep all saved 1-mm/h object identities."""
import argparse
from collections import Counter
from datetime import datetime, timedelta, timezone
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from audit_reference_stream import inspect
from audit_seed_strength import verify_manifest_coverage
from compare_frozen_extents import digest_frame, manifest
from compare_rain_thresholds import dump, reference_equivalence, sha_file
from frame_progress import FrameProgress
from run_npy_validation import disk, write_graph
from step.envelope import expand_seeded_envelope
from step.identification import identify
from step.tracking import load_tracking_state, save_tracking_state, track_with_graph
from validate_extent_holdouts import deadline, stop_group

GEOMETRIES = {'seed_control': None, 'extent_0p5': .5, 'extent_0p1': .1}
NOTE = ('Catalog entry remains the saved 1-mm/h seed definition; all frame/local-label '
        'identities are retained. Only experimental tracking geometry/statistics change. '
        'The expanded support jointly changes overlap, centroid, area and mean intensity; '
        'this is not an attribution to any single feature. Original graph/catalog files '
        'are immutable and not replaced. More links or longer branches do not prove '
        'accuracy. Expanded branch numbers are independent of baseline branch numbers.')


def compare_links(reference, current):
    before, old_objects, old_edges, _, _ = inspect(reference)
    after, new_objects, new_edges, _, _ = inspect(current)
    def identities(objects):
        return {int(r['node_id']): (int(r['time']), int(r['local_label'])) for r in objects}
    old_ids, new_ids = identities(old_objects), identities(new_objects)
    if old_ids != new_ids:
        raise AssertionError('Fixed seed node/frame/local-label identities changed')
    def keyed(edges, ids):
        return {ids[int(r['parent_node_id'])]+ids[int(r['child_node_id'])]: r for r in edges}
    old, new = keyed(old_edges, old_ids), keyed(new_edges, new_ids)
    added = [dict(pair=list(k), event=new[k]['event'], score=float(new[k]['score']))
             for k in sorted(new.keys()-old.keys())]
    removed = [dict(pair=list(k), event=old[k]['event'], score=float(old[k]['score']))
               for k in sorted(old.keys()-new.keys())]
    changed = [dict(pair=list(k), before=old[k]['event'], after=new[k]['event'])
               for k in sorted(old.keys() & new.keys()) if old[k]['event'] != new[k]['event']]
    deltas = [abs(float(old[k]['score'])-float(new[k]['score'])) for k in old.keys() & new.keys()]
    return dict(fixed_node_identities_equal=True, before=before, after=after,
        retained_links=len(old.keys() & new.keys()), added=len(added), removed=len(removed),
        changed_event=len(changed), added_edge_types=dict(Counter(r['event'] for r in added)),
        removed_edge_types=dict(Counter(r['event'] for r in removed)),
        retained_link_score_delta_median=float(np.median(deltas)) if deltas else None,
        note=NOTE), dict(added=added, removed=removed, changed_event=changed)


def run(extents, output, chunk_frames=6, workers=4, finish_before=None):
    if output.exists():
        raise FileExistsError(output)
    if chunk_frames < 1 or workers < 1:
        raise ValueError('Positive workers and chunk size required')
    if not (extents / 'SUCCESS').is_file():
        raise ValueError('Successful fixed-seed extent run required')
    started = time.perf_counter()
    config = json.loads((extents / 'configuration.json').read_text())
    source, metadata, reference = (Path(config[k]) for k in ('prepared_input', 'metadata', 'reference'))
    if sha_file(source) != config['input_sha256'] or sha_file(metadata) != config['metadata_sha256']:
        raise ValueError('Prepared input changed')
    before = manifest(reference)
    added_coverage = verify_manifest_coverage(config['source_graph_files_sha256'], before, reference)
    extent_manifest = {name: sha_file(extents / name) for name in
        ('SUCCESS', 'configuration.json', 'summary.json', 'seed_frame_hashes.json')}
    hashes = json.loads((extents / 'seed_frame_hashes.json').read_text())
    rain = np.load(source, mmap_mode='r')
    meta = json.loads(metadata.read_text())
    structure, objects, _, _, _ = inspect(reference)
    if len(rain) != len(hashes) or len(rain) != structure['frames']:
        raise ValueError('Saved sequence length differs')
    nodes = {}
    for row in objects:
        nodes.setdefault(int(row['time']), {})[int(row['local_label'])] = row
    output.mkdir(parents=True)
    dump(output / 'PLAN.json', dict(registered_at=datetime.now(timezone.utc).isoformat(),
        geometries=GEOMETRIES, preset=config['preset'], source=str(extents.resolve()),
        chunk_frames=chunk_frames, workers=workers, input_sha256=config['input_sha256'],
        metadata_sha256=config['metadata_sha256'], source_graph_files_sha256=before,
        extent_source_sha256=extent_manifest, manifest_coverage_added=added_coverage,
        code_sha256=sha_file(Path(__file__)), note=NOTE))
    states = {name: None for name in GEOMETRIES}
    for name in GEOMETRIES:
        (output / name).mkdir()
    footprint = disk(config['preset']['identification']['bridge_radius_cells'])
    kwargs = dict(config['preset']['tracking'])
    kwargs['km'] = kwargs.pop('max_displacement')
    timings = {name: 0. for name in GEOMETRIES}
    progress = FrameProgress(len(rain))
    progress.unit, progress.verb = 'hours', 'compared'
    for offset in range(0, len(rain), chunk_frames):
        if finish_before and datetime.now(timezone.utc) >= finish_before-timedelta(minutes=10):
            raise TimeoutError('Deadline reserve reached; incomplete experiment not successful')
        data = np.array(rain[offset:offset+chunk_frames], copy=True)
        cores = identify(data, footprint, workers=workers, threshold=1,
                         min_size=config['preset']['identification']['min_size_cells'])
        for i, seeds in enumerate(cores):
            if digest_frame(seeds) != hashes[offset+i]:
                raise AssertionError('Frozen seed hash changed')
            counts = np.bincount(seeds.ravel())
            if (set(np.flatnonzero(counts[1:])+1) != set(nodes[offset+i]) or
                any(counts[label] != int(row['area_cells']) for label, row in nodes[offset+i].items())):
                raise AssertionError('Frozen seed labels/areas changed')
        for name, threshold in GEOMETRIES.items():
            progress.update(offset, phase=f'{name}: geometry and tracking; original graph untouched')
            t = time.perf_counter()
            geometry = cores if threshold is None else np.stack([
                expand_seeded_envelope(frame, seeds, threshold) for frame, seeds in zip(data, cores)])
            for seeds, labels in zip(cores, geometry):
                if not np.array_equal(labels[seeds > 0], seeds[seeds > 0]):
                    raise AssertionError('Tracking geometry changed seed labels')
                if set(np.unique(labels[labels > 0])) != set(np.unique(seeds[seeds > 0])):
                    raise AssertionError('Tracking geometry changed object entry')
            tracked, graph, states[name] = track_with_graph(geometry, data,
                state=states[name], return_state=True, start_time=offset,
                sequence_id='geometry_'+config['input_sha256'][:16]+'_'+name, **kwargs)
            part = output / name / f'chunk_{offset:06d}'
            part.mkdir()
            write_graph(graph, part)
            save_tracking_state(states[name], part / 'state.json')
            states[name] = load_tracking_state(part / 'state.json')
            dump(part / 'metadata.json', dict(timestamps=meta['timestamps'][offset:offset+len(data)],
                dt_hours=1., normalized_units='mm/h', grid_attributes=meta['grid'],
                crop_origin_yx=[0, 0], tracking_geometry=name))
            dump(part / 'frame_hashes.json', dict(
                identified_frame_sha256=[digest_frame(frame) for frame in cores],
                tracking_geometry_frame_sha256=[digest_frame(frame) for frame in geometry],
                tracked_frame_sha256=[digest_frame(frame) for frame in tracked]))
            timings[name] += time.perf_counter()-t
            del graph, tracked, geometry
        progress.update(offset+len(data), phase='all three tracking geometries completed for this chunk')
    for name in GEOMETRIES:
        dump(output / name / 'SUCCESS', dict(scope='Experimental tracking geometry, not promoted'))
    # Control must reproduce the complete saved graph, not merely its node count.
    control = reference_equivalence(output / 'seed_control', reference)
    comparisons = {}
    for name in ('extent_0p5', 'extent_0p1'):
        comparisons[name], changes = compare_links(reference, output / name)
        dump(output / (name+'_edge_changes.json'), changes)
    if (manifest(reference) != before or
        {name: sha_file(extents / name) for name in extent_manifest} != extent_manifest or
        sha_file(source) != config['input_sha256'] or sha_file(metadata) != config['metadata_sha256']):
        raise AssertionError('Frozen sources changed')
    summary = dict(frames=len(rain), first_interval_end=meta['timestamps'][0],
        baseline_control=control, comparisons=comparisons, stage_wall_seconds=timings,
        wall_seconds=time.perf_counter()-started, source_files_unchanged=True, note=NOTE)
    dump(output / 'summary.json', summary)
    dump(output / 'SUCCESS', dict(scope='Controlled sensitivity, not tracking accuracy'))
    progress.update(phase='SUCCESS: experimental graphs saved; no default changed')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('extents', type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--chunk-frames', type=int, default=6)
    parser.add_argument('--finish-before', required=True, type=deadline)
    parser.add_argument('--timeout-seconds', type=float, default=1200)
    parser.add_argument('--internal-worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.internal_worker:
        run(args.extents, args.output_dir, args.chunk_frames, args.workers, args.finish_before)
        return
    remaining = (args.finish_before-timedelta(minutes=10)-datetime.now(timezone.utc)).total_seconds()
    if remaining <= 0 or args.timeout_seconds <= 0:
        raise ValueError('No time remains before reporting reserve / invalid timeout')
    process = subprocess.Popen([sys.executable, '-u', str(Path(__file__).resolve()),
                               *sys.argv[1:], '--internal-worker'], start_new_session=True)
    try:
        code = process.wait(timeout=min(remaining, args.timeout_seconds))
    except BaseException:
        stop_group(process)
        raise
    if code:
        raise SystemExit(code)


if __name__ == '__main__':
    main()
