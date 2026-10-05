"""Retrospective strong-phase sensitivity on an immutable 1-mm/h lineage graph."""
import argparse
from collections import defaultdict
from datetime import datetime, timedelta, timezone
import csv
import json
from pathlib import Path
import time

import numpy as np
from scipy import ndimage

from audit_reference_stream import inspect
from compare_frozen_extents import digest_frame, manifest
from compare_rain_thresholds import dump, sha_file, write_rows
from frame_progress import FrameProgress
from run_npy_validation import disk
from step.identification import identify_frame
from summarize_extent_holdouts import CHECKS
from validate_extent_holdouts import deadline

THRESHOLDS = (1., 2., 5.)
AREAS = (1, 4, 16)
DURATIONS = (1, 3)
NOTE = ('Descriptive strong-rain phase audit, not new identification/tracking or convection '
        'classification. Area means largest contiguous 8-neighbor raw strong patch inside '
        'a saved 1-mm/h object, not dilated area or total scattered strong pixels. '
        'Branch phase is retrospective over observed event-ended branches; split/merge '
        'boundaries truncate branches and gaps interrupt consecutive qualification. '
        'Window endpoints censor lifetimes. Weak-only seed-free phases remain absent. '
        'Selection does not redistribute excluded watershed pixels to surviving seeds. '
        'Rain sums use equal grid weights and supplied RAINNC only.')


def verify_manifest_coverage(recorded, current, reference):
    """Keep every old signature; allow only two explicitly known coverage additions."""
    changed = [name for name, digest in recorded.items() if current.get(name) != digest]
    added = set(current)-set(recorded)
    allowed = {str((reference.parent / name).resolve()) for name in
               ('frame_accounting.csv', 'c_envelope_objects.csv')}
    if changed or added-allowed:
        raise ValueError('Frozen graph source changed or has unexpected unsigned files')
    return sorted(added)


def strong_patches(rain, seeds, threshold):
    if (rain.shape != seeds.shape or rain.ndim != 2 or not np.isfinite(threshold) or threshold < 1 or
        not np.issubdtype(seeds.dtype, np.integer) or np.any(seeds < 0)):
        raise ValueError('Requires matched 2D seeds/rain and threshold >=1')
    seed_mask = seeds > 0
    if np.any(seed_mask & (~np.isfinite(rain) | (rain < 1))):
        raise ValueError('Saved seed support must be finite and >=1 mm/h')
    strong = seed_mask & (rain >= threshold)
    size = int(seeds.max())+1
    counts = np.bincount(seeds[strong], minlength=size)
    components, count = ndimage.label(strong, np.ones((3, 3), bool))
    maximum = np.zeros(size, dtype=np.int64)
    patch_count = np.zeros(size, dtype=np.int64)
    if count:
        active = np.flatnonzero(components)
        ids = components.ravel()[active]
        first = np.full(count+1, rain.size, dtype=np.int64)
        np.minimum.at(first, ids, active)
        owners = np.zeros(count+1, dtype=seeds.dtype)
        owners[1:] = seeds.ravel()[first[1:]]
        if np.any(seeds.ravel()[active] != owners[ids]):
            raise AssertionError('Contiguous strong patch crosses saved object IDs')
        sizes = np.bincount(ids, minlength=count+1)
        np.maximum.at(maximum, owners[1:], sizes[1:])
        patch_count = np.bincount(owners[1:], minlength=size)
    return counts, maximum, patch_count


def longest_run(values):
    """Values are (absolute hourly frame, qualifies); missing hours break a run."""
    best = current = 0
    previous = None
    for stamp, qualifies in sorted(values):
        current = (current+1 if previous is not None and stamp == previous+1 else 1) if qualifies else 0
        best = max(best, current)
        previous = stamp
    return best


def policies(records, masses, all_rain, assigned):
    result = []
    branches = defaultdict(list)
    nodes = {}
    for row in records:
        key = (row['threshold_mm_h'], row['branch_id'])
        branches[key].append(row)
        nodes[row['node_id']] = row['branch_id']
    for threshold in THRESHOLDS:
        rows = [r for r in records if r['threshold_mm_h'] == threshold]
        for area in AREAS:
            snapshot = {r['node_id'] for r in rows if r['largest_patch_cells'] >= area}
            runs = {branch: longest_run([(r['frame'], r['largest_patch_cells'] >= area)
                                       for r in values])
                    for (t, branch), values in branches.items() if t == threshold}
            selections = [('snapshot', None, snapshot)]
            for duration in DURATIONS:
                qualified = {branch for branch, length in runs.items() if length >= duration}
                selections.append(('branch_phase', duration,
                    {r['node_id'] for r in rows if r['branch_id'] in qualified}))
            for policy, duration, selected in selections:
                rain = {str(t): sum(masses[t][node] for node in selected) for t in (.5, .1)}
                result.append(dict(threshold_mm_h=threshold, min_patch_cells=area,
                    selection=policy, min_consecutive_hours=duration,
                    qualifying_snapshots=len(snapshot), retained_snapshots=len(selected),
                    retained_branches=len({nodes[node] for node in selected}),
                    retained_extent_rain=rain,
                    retained_fraction_of_all_rain={t: value/all_rain if all_rain else None
                                                  for t, value in rain.items()},
                    retained_fraction_of_assigned_rain={t: value/assigned[t] if assigned[t] else None
                                                       for t, value in rain.items()}))
    return result


def branch_duration_context(objects, edges, masses, frames):
    """Describe overlapping censoring flags, not reasons for physical short life."""
    groups = defaultdict(list)
    by_node = {int(row['node_id']): row for row in objects}
    for row in objects:
        groups[int(row['branch_id'])].append((int(row['time']), int(row['node_id'])))
    event_branches, gap_branches = set(), set()
    for row in edges:
        endpoints = (int(row['parent_node_id']), int(row['child_node_id']))
        if row['event'] in ('split', 'merge', 'complex'):
            event_branches.update(int(by_node[n]['branch_id']) for n in endpoints)
        if row['event'] == 'gap_continue':
            gap_branches.update(int(by_node[n]['branch_id']) for n in endpoints)
    rejected = {b for b, values in groups.items() if longest_run([(t, True) for t, _ in values]) < 3}
    boundary = {b for b, values in groups.items() if any(t in (0, frames-1) for t, _ in values)}
    flagged = boundary | event_branches | gap_branches
    subsets = dict(all_short_consecutive_branches=rejected,
        touching_window_boundary=rejected & boundary,
        touching_event_endpoint=rejected & event_branches,
        containing_gap=rejected & gap_branches,
        without_these_flags=rejected-flagged)
    total = sum(masses.values())
    return dict(min_consecutive_observed_hours=3,
        subsets={name: dict(branches=len(selected),
            object_snapshots=sum(len(groups[b]) for b in selected),
            assigned_0p1_rain_fraction=sum(masses[node] for b in selected for _, node in groups[b])/total
                if total else None) for name, selected in subsets.items()},
        note='Flags overlap and must not be added. Window contact is possible censoring, '
             'not confirmed physical birth/decay. Event endpoints are algorithmic branch '
             'boundaries. Unflagged rejection is not proof of a short physical storm.')


def run(extents, output, finish_before=None):
    if output.exists():
        raise FileExistsError(output)
    started = time.perf_counter()
    if not (extents / 'SUCCESS').is_file():
        raise ValueError('Successful fixed-seed extent output required')
    summary = json.loads((extents / 'summary.json').read_text())
    if not all(summary.get(key) is True for key in CHECKS):
        raise ValueError('Extent source integrity checks did not pass')
    if any(summary['extents'][t]['assigned_rain'] <= 0 for t in ('0.5', '0.1')):
        raise ValueError('Rain-fraction audit requires positive assigned rain at both extents')
    config = json.loads((extents / 'configuration.json').read_text())
    source, metadata, reference = (Path(config[k]) for k in ('prepared_input', 'metadata', 'reference'))
    if sha_file(source) != config['input_sha256'] or sha_file(metadata) != config['metadata_sha256']:
        raise ValueError('Prepared input checksums changed')
    before = manifest(reference)
    coverage_added = verify_manifest_coverage(config['source_graph_files_sha256'], before, reference)
    source_files = ('SUCCESS', 'summary.json', 'configuration.json', 'extent_objects.csv',
                    'seed_frame_hashes.json')
    before_extent = {name: sha_file(extents / name) for name in source_files}
    structure, objects, edges, _, _ = inspect(reference)
    nodes = defaultdict(dict)
    for row in objects:
        nodes[int(row['time'])][int(row['local_label'])] = row
    masses = {.5: {}, .1: {}}
    with (extents / 'extent_objects.csv').open(newline='') as handle:
        for row in csv.DictReader(handle):
            t, nid = float(row['threshold_mm_h']), int(row['core_node_id'])
            if nid in masses[t]:
                raise ValueError('Duplicate extent object')
            masses[t][nid] = float(row['envelope_rain'])
    node_ids = {int(row['node_id']) for row in objects}
    assigned = {t: summary['extents'][t]['assigned_rain'] for t in ('0.5', '0.1')}
    for threshold, values in masses.items():
        if set(values) != node_ids or not np.isclose(sum(values.values()), assigned[str(threshold)], rtol=1e-12):
            raise ValueError('Extent/node accounting differs')
    rain = np.load(source, mmap_mode='r')
    hashes = json.loads((extents / 'seed_frame_hashes.json').read_text())
    meta = json.loads(metadata.read_text())
    if len(rain) != len(hashes) or len(rain) != structure['frames']:
        raise ValueError('Frame extent mismatch')
    output.mkdir(parents=True)
    dump(output / 'PLAN.json', dict(registered_at=datetime.now(timezone.utc).isoformat(),
        thresholds_mm_h=THRESHOLDS, min_patch_cells=AREAS, consecutive_hours=DURATIONS,
        source=str(extents.resolve()), source_artifacts_sha256=before_extent,
        source_graph_files_sha256=before, finish_before=finish_before.isoformat() if finish_before else None,
        manifest_coverage_added_since_extent_run=coverage_added,
        code_sha256=sha_file(Path(__file__)), note=NOTE))
    progress = FrameProgress(len(rain))
    progress.unit, progress.verb = 'hours', 'audited'
    footprint = disk(config['preset']['identification']['bridge_radius_cells'])
    records = []
    for frame in range(len(rain)):
        if finish_before and datetime.now(timezone.utc) >= finish_before-timedelta(minutes=10):
            raise TimeoutError('Reporting reserve reached; incomplete audit not marked successful')
        current = np.array(rain[frame], copy=True)
        seeds = identify_frame(current, footprint, threshold=1,
            min_size=config['preset']['identification']['min_size_cells'])
        if digest_frame(seeds) != hashes[frame]:
            raise AssertionError('Frozen seed hash differs')
        counts = np.bincount(seeds.ravel())
        if set(np.flatnonzero(counts[1:])+1) != set(nodes[frame]):
            raise AssertionError('Frame seed labels differ from catalog')
        previous_strong = previous_largest = None
        for threshold in THRESHOLDS:
            strong, largest, patches = strong_patches(current, seeds, threshold)
            if previous_strong is not None and (np.any(strong > previous_strong) or
                                               np.any(largest > previous_largest)):
                raise AssertionError('Stronger threshold increased raw support or largest patch')
            previous_strong, previous_largest = strong, largest
            for label, node in nodes[frame].items():
                if counts[label] != int(node['area_cells']):
                    raise AssertionError('Frozen seed area differs')
                records.append(dict(frame=frame, timestamp=meta['timestamps'][frame],
                    node_id=int(node['node_id']), branch_id=int(node['branch_id']),
                    local_label=label, threshold_mm_h=threshold, seed_cells=int(counts[label]),
                    strong_cells=int(strong[label]), largest_patch_cells=int(largest[label]),
                    strong_patch_count=int(patches[label]), max_rate_mm_h=float(node['max_intensity'])))
        progress.update(frame+1, phase='strength/area statistics; graph unchanged')
    results = policies(records, masses, summary['extents']['0.1']['all_provided_rain'], assigned)
    control = next(r for r in results if r['selection'] == 'snapshot' and
                   r['threshold_mm_h'] == 1 and r['min_patch_cells'] == 1)
    if control['retained_snapshots'] != len(objects) or any(
        not np.isclose(v, 1, rtol=1e-12) for v in control['retained_fraction_of_assigned_rain'].values()):
        raise AssertionError('1-mm/h, one-pixel control does not preserve all assigned rain')
    if (manifest(reference) != before or
        {name: sha_file(extents / name) for name in source_files} != before_extent or
        sha_file(source) != config['input_sha256'] or sha_file(metadata) != config['metadata_sha256']):
        raise AssertionError('Immutable sources changed during audit')
    context = branch_duration_context(objects, edges, masses[.1], len(rain))
    duration_control = next(r for r in results if r['selection'] == 'branch_phase' and
        r['threshold_mm_h'] == 1 and r['min_patch_cells'] == 1 and r['min_consecutive_hours'] == 3)
    rejected_fraction = context['subsets']['all_short_consecutive_branches']['assigned_0p1_rain_fraction']
    if node_ids and not np.isclose(rejected_fraction+
        duration_control['retained_fraction_of_assigned_rain']['0.1'], 1, rtol=1e-12):
        raise AssertionError('Duration-context rain does not conserve control selection')
    report = dict(frames=len(rain), first_interval_end=meta['timestamps'][0],
        last_interval_end=meta['timestamps'][-1], structure=structure, policies=results,
        duration_context=context, manifest_coverage_added_since_extent_run=coverage_added,
        exact_seed_hashes_equal=True, graph_and_extent_sources_unchanged=True,
        baseline_preserves_all_assigned_rain=True, wall_seconds=time.perf_counter()-started, note=NOTE)
    write_rows(output / 'object_strength.csv', records)
    dump(output / 'summary.json', report)
    lines = ['# Frozen-graph seed-strength audit', '', NOTE, '',
        '| Rate mm/h | Contiguous patch cells | Selection | Consecutive hours | Retained snapshots | Retained branches | 0.5 / 0.1 extent rain retained (% of assigned) |',
        '|---|---:|---|---:|---:|---:|---:|']
    for row in results:
        fractions = row['retained_fraction_of_assigned_rain']
        lines.append(f"| {row['threshold_mm_h']} | {row['min_patch_cells']} | {row['selection']} | "
            f"{row['min_consecutive_hours'] or '--'} | {row['retained_snapshots']} | "
            f"{row['retained_branches']} | {100*fractions['0.5']:.3f} / {100*fractions['0.1']:.3f} |")
    context = report['duration_context']
    lines += ['', '## Three-hour observed-branch control', '', context['note'], '',
        '| Rejected subset | Branches | Snapshots | % assigned 0.1 rain |',
        '|---|---:|---:|---:|']
    for name, values in context['subsets'].items():
        lines.append(f"| {name} | {values['branches']} | {values['object_snapshots']} | "
                     f"{100*values['assigned_0p1_rain_fraction']:.3f} |")
    (output / 'REPORT.md').write_text('\n'.join(lines)+'\n')
    dump(output / 'SUCCESS', dict(scope='Integrity and sensitivity audit, not scientific certification'))
    progress.update(phase='SUCCESS: audit recorded; defaults unchanged')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('extents', type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--finish-before', type=deadline,
                        help='Timezone-aware cutoff; reserve final ten minutes for reporting')
    args = parser.parse_args()
    run(args.extents, args.output_dir, args.finish_before)
