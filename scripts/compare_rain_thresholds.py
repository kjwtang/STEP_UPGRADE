#!/usr/bin/env python3
"""Paired rain-threshold sensitivity; fixed tracking, not storm accuracy validation.

A tracks >=1 mm/h. B tracks >=0.5 mm/h with the identical tracking settings.
C partitions connected >=0.5 rain among A seeds and reuses A's graph: it does
not track envelope centroids. Whole weak components are a bridge diagnostic,
not an additional tracking variant. Inputs are already normalized hourly rates.
"""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from step.envelope import identify_envelope
from step.identification import identify
from step.tracking import track_with_graph, save_tracking_state, load_tracking_state
from audit_reference_stream import inspect, frame_hashes, input_identity
from frame_progress import FrameProgress
from run_npy_validation import disk, write_graph

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PRESET = ROOT / 'configs/rainfall_lineage_4km_1h_reference_v1.json'
LIMITATIONS = ('Sensitivity, not accuracy. A/B labels and IDs are not directly comparable. '
    'Only bijective shared-core associations enter projected edge comparisons; multicore '
    'groupings are reported separately, not called false merges. C keeps A tracking and '
    'does not validate tracking expanded masks. Operational rain seeds are not diagnosed '
    'convective cores. Rain sums use equal grid cells and provided RAINNC only. '
    'One 1996 window is not independent historical/future validation. Timings are single '
    'warm-input measurements, not a controlled speed comparison.')


def sha_file(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False))


def write_rows(path, rows):
    if not rows:
        dump(path.with_suffix('.json'), [])
        return
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def core_memberships(core, lower, frame):
    """Sparse shared-pixel mapping; never equate labels across definitions."""
    if np.any((core > 0) & (lower == 0)):
        raise AssertionError('Lower threshold lost a retained reference pixel')
    mask = core > 0
    base = int(lower.max()) + 1
    codes, counts = np.unique(core[mask].astype(np.int64) * base + lower[mask],
                              return_counts=True)
    return [dict(frame=frame, core_label=int(code // base), lower_label=int(code % base),
                 shared_core_pixels=int(count)) for code, count in zip(codes, counts)]


def accounting(rain, core, lower, partition, system):
    valid = np.isfinite(rain)
    if np.any(rain[valid] < 0):
        raise ValueError('Negative hourly rate; preprocess cumulative resets explicitly')
    if not np.array_equal(partition[core > 0], core[core > 0]):
        raise AssertionError('Partition changed core ownership')
    if not np.array_equal(partition > 0, system > 0):
        raise AssertionError('Partition/system assigned support differs')
    if np.any((partition > 0) & (~valid | (rain < .5))):
        raise AssertionError('Envelope crossed a dry/missing barrier')
    def total(mask):
        return float(np.sum(rain[mask], dtype=np.float64))
    assigned = partition > 0
    unassigned = (lower > 0) & ~assigned
    result = dict(valid_cells=int(valid.sum()), missing_cells=int((~valid).sum()),
        all_provided_rain=total(valid), a_cells=int((core > 0).sum()),
        b_cells=int((lower > 0).sum()), c_cells=int(assigned.sum()),
        a_rain=total(core > 0), b_rain=total(lower > 0), c_rain=total(assigned),
        unassigned_lower_cells=int(unassigned.sum()), unassigned_lower_rain=total(unassigned))
    if not np.isclose(result['c_rain'] + result['unassigned_lower_rain'],
                      result['b_rain'], rtol=1e-12, atol=1e-8):
        raise AssertionError('Seeded/unassigned rain accounting does not conserve B rain')
    return result


def envelope_geometry(rain, core, partition, nodes, frame):
    """Intensity-weighted centroid drift; C never uses these centroids to link."""
    active = np.flatnonzero(partition)
    labels = partition.ravel()[active]
    values = rain.ravel()[active].astype(np.float64)
    n = int(core.max()) + 1
    counts = np.bincount(labels, minlength=n)
    mass = np.bincount(labels, weights=values, minlength=n)
    sy = np.bincount(labels, weights=values * (active // rain.shape[1]), minlength=n)
    sx = np.bincount(labels, weights=values * (active % rain.shape[1]), minlength=n)
    rows = []
    for label, node in sorted(nodes.items()):
        y, x = sy[label] / mass[label], sx[label] / mass[label]
        rows.append(dict(frame=frame, core_label=label, core_node_id=node.node_id,
            branch_id=node.branch_id, core_cells=node.object.area,
            envelope_cells=int(counts[label]), envelope_rain=float(mass[label]),
            envelope_to_core_cell_ratio=float(counts[label] / node.object.area),
            centroid_drift_cells=float(np.hypot(y-node.object.centroid[0], x-node.object.centroid[1]))))
    return rows


def projected_edges(a_objects, a_edges, b_objects, b_edges, memberships):
    """Only one-to-one shared-core endpoint correspondences; no Cartesian truth."""
    by_a, by_b = defaultdict(set), defaultdict(set)
    for row in memberships:
        a, b = (row['frame'], row['core_label']), (row['frame'], row['lower_label'])
        by_a[a].add(b)
        by_b[b].add(a)
    mapping = {b: next(iter(values)) for b, values in by_b.items()
               if len(values) == 1 and len(by_a[next(iter(values))]) == 1}
    comparable = set(mapping.values())
    def edges(rows, objects, translate=None):
        nodes = {int(r['node_id']): (int(r['time']), int(r['local_label'])) for r in objects}
        result = {}
        excluded = 0
        for row in rows:
            p, c = nodes[int(row['parent_node_id'])], nodes[int(row['child_node_id'])]
            if translate is not None:
                if p not in translate or c not in translate:
                    excluded += 1
                    continue
                p, c = translate[p], translate[c]
            elif p not in comparable or c not in comparable:
                excluded += 1
                continue
            result[p + c] = row['event']
        return result, excluded
    old, a_excluded = edges(a_edges, a_objects)
    new, b_excluded = edges(b_edges, b_objects, mapping)
    return dict(bijectively_mapped_core_snapshots=len(mapping),
        excluded_a_edges=a_excluded, excluded_b_edges=b_excluded,
        comparable_a_edges=len(old), comparable_b_edges=len(new),
        retained=len(old.keys() & new.keys()), added=len(new.keys() - old.keys()),
        removed=len(old.keys() - new.keys()),
        changed_event=sum(old[k] != new[k] for k in old.keys() & new.keys()),
        note='Shared >=1-mm/h pixel association, not accuracy or whole-graph equivalence.'), dict(
            added=[dict(pair=list(k), event=new[k]) for k in sorted(new.keys() - old.keys())],
            removed=[dict(pair=list(k), event=old[k]) for k in sorted(old.keys() - new.keys())],
            changed_event=[dict(pair=list(k), before=old[k], after=new[k])
                           for k in sorted(old.keys() & new.keys()) if old[k] != new[k]])


def catalog_metrics(objects, edges, frames):
    branches = defaultdict(list)
    incoming = {int(r['child_node_id']) for r in edges}
    total = sum(float(r['precipitation_sum']) for r in objects)
    small = [r for r in objects if int(r['area_cells']) < 16]
    eligible = [r for r in objects if int(r['time']) > 0]
    linked = [r for r in eligible if int(r['node_id']) in incoming]
    for row in objects:
        branches[row['branch_id']].append(int(row['time']))
    spans = [max(values)-min(values)+1 for values in branches.values()]
    denominator = sum(float(r['precipitation_sum']) for r in eligible)
    return dict(small_object_snapshots=len(small),
        small_object_fraction=len(small)/len(objects) if objects else None,
        small_object_labeled_rain_fraction=sum(float(r['precipitation_sum']) for r in small)/total if total else None,
        one_snapshot_branch_fraction=sum(len(v) == 1 for v in branches.values())/len(branches) if branches else None,
        branch_span_hours_median=float(np.median(spans)) if spans else None,
        branch_span_hours_max=max(spans, default=None),
        branches_touching_time_boundary=sum(0 in v or frames-1 in v for v in branches.values()),
        incoming_object_association_fraction=len(linked)/len(eligible) if eligible else None,
        incoming_labeled_rain_association_fraction=sum(float(r['precipitation_sum']) for r in linked)/denominator if denominator else None,
        note='Observed branch spans include gaps and event cuts; association coverage is not accuracy.')


def reference_equivalence(source, reference):
    """Require matching saved chronology/grid as well as every scientific output."""
    _, o, e, v, r = inspect(reference)
    _, ao, ae, av, ar = inspect(source)
    identity, other_identity = input_identity(source), input_identity(reference)
    # The preparation metadata includes shape among attributes; old stream
    # metadata stores shape separately. Verify that redundant shape, then
    # compare the identical physical attributes and checkpoint grid shape.
    for value in (identity, other_identity):
        attributes = value['grid']['grid_attributes']
        shape = attributes.pop('shape', None) if attributes is not None else None
        if shape is not None and shape != value['grid']['processed_grid_shape']:
            raise ValueError('Saved attribute shape differs from checkpoint grid')
    def hashes(root, key):
        parts = sorted(root.glob('chunk_*'))
        if all((p / 'frame_hashes.json').exists() for p in parts):
            return [h for p in parts for h in json.loads((p / 'frame_hashes.json').read_text())[key]]
        return frame_hashes(root, key)
    equality = dict(objects_equal=ao == o, full_edges_scores_equal=ae == e,
        events_equal=av == v, final_roots_equal=ar == r,
        saved_timestamps_equal=identity['timestamps'] is not None and identity['timestamps'] == other_identity['timestamps'],
        saved_grid_units_cadence_equal=identity['grid'] == other_identity['grid'],
        identified_raster_hashes_equal=hashes(source, 'identified_frame_sha256') == hashes(reference, 'identified_frame_sha256'),
        tracked_raster_hashes_equal=hashes(source, 'tracked_frame_sha256') == hashes(reference, 'tracked_frame_sha256'))
    if not all(equality.values()):
        raise AssertionError(f'A frozen-reference reproduction failed: {equality}')
    return equality


def render(frame, stamp, rain, a, b, partition, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.8), layout='constrained')
    im = axes[0].imshow(rain, origin='lower', vmin=0, vmax=10, cmap='Blues')
    fig.colorbar(im, ax=axes[0], label='mm/h')
    axes[0].set_title('Rain; display capped at 10')
    for ax, labels, title in zip(axes[1:], (a, b, partition),
            ('A: 1.0 tracked branches', 'B: 0.5 tracked branches', 'C: A branches + 0.5 envelope')):
        ax.imshow(np.ma.masked_where(labels == 0, labels % 20), origin='lower',
                  cmap='tab20', vmin=0, vmax=19, interpolation='nearest')
        ax.set_title(title)
    for ax in axes:
        ax.set(xlabel='x (full-grid cells)', ylabel='y (full-grid cells)')
    fig.suptitle(f'{stamp} | frame {frame}; A/B IDs independent; colors repeat')
    fig.savefig(output / f'frame_{frame:03d}.png', dpi=140)
    plt.close(fig)


def run(source, metadata_path, output, preset_path=DEFAULT_PRESET, workers=4,
        chunk_frames=6, reference=None, plots=True):
    if workers < 1 or chunk_frames < 1:
        raise ValueError('workers and chunk-frames must be positive')
    if output.exists():
        raise FileExistsError(output)
    preset = json.loads(preset_path.read_text())
    if preset['identification']['threshold_mm_h'] != 1 or preset['identification']['min_size_cells'] != 1:
        raise ValueError('This paired experiment requires 1-mm/h reference and min-size=1')
    if preset['input_contract'] != dict(grid_km=4, dt_hours=1, rain_units_after_normalization='mm/h'):
        raise ValueError('Experiment requires the explicit 4-km/1-hour rate contract')
    rain = np.load(source, mmap_mode='r')
    meta = json.loads(metadata_path.read_text())
    stamps = np.asarray(meta['timestamps'], dtype='datetime64[s]')
    if (rain.ndim != 3 or len(rain) < 2 or not np.issubdtype(rain.dtype, np.floating)
            or len(stamps) != len(rain) or np.isnat(stamps).any()
            or not np.all(np.diff(stamps) == np.timedelta64(1, 'h'))):
        raise ValueError('Require floating (time,y,x) rates and complete hourly timestamps')
    grid = meta['grid']
    if (list(rain.shape[1:]) != grid['shape'] or grid['DX'] != 4000 or grid['DY'] != 4000):
        raise ValueError('Prepared rate shape / 4-km grid mismatch')
    source_sha = sha_file(source)
    output.mkdir(parents=True)
    modes = dict(A_1mm=1., B_0p5mm=.5)
    settings = dict(input_sha256=source_sha, metadata_sha256=sha_file(metadata_path),
        units='mm/h', dt_hours=1, workers=workers, chunk_frames=chunk_frames,
        prepared_input=str(source.resolve()), preset=preset, limitations=LIMITATIONS,
        variants=dict(A_identification_threshold_mm_h=1., B_identification_threshold_mm_h=.5,
                      C_seed_threshold_mm_h=1., C_extent_threshold_mm_h=.5,
                      C_tracking_graph='Reused A, not newly tracked'),
        python_version=sys.version, dependency_versions={name: version(name)
            for name in ('numpy', 'scipy', 'scikit-image')},
        source_sha256={p: sha_file(ROOT / p) for p in ('step/tracking.py',
            'step/identification.py', 'step/envelope.py', 'scripts/compare_rain_thresholds.py')},
        git_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip())
    dump(output / 'configuration.json', settings)
    states = {name: None for name in modes}
    times = {name: dict(identification_seconds=0., tracking_seconds=0., output_seconds=0.) for name in modes}
    for name, threshold in modes.items():
        (output / name).mkdir()
        dump(output / name / 'settings.json', dict(threshold=threshold, tracking=preset['tracking']))
    kwargs = dict(preset['tracking'])
    kwargs['km'] = kwargs.pop('max_displacement')
    sequence = 'paired_threshold_' + source_sha[:16]
    footprint = disk(preset['identification']['bridge_radius_cells'])
    frame_rows, mappings, geometries, bridge_rows = [], [], [], []
    c_seconds = 0.
    started = time.perf_counter()
    progress = FrameProgress(len(rain))
    selected = set(np.linspace(0, len(rain)-1, min(4, len(rain)), dtype=int).tolist())
    succeeded = False
    try:
        for offset in range(0, len(rain), chunk_frames):
            data = np.array(rain[offset:offset+chunk_frames], copy=True)
            if np.any(data[np.isfinite(data)] < 0):
                raise ValueError('Negative hourly rate; do not guess cumulative reset handling')
            labels, tracked, graphs = {}, {}, {}
            for name, threshold in modes.items():
                progress.update(offset, phase=f'{name}: identification / tracking')
                t = time.perf_counter()
                labels[name] = identify(data, footprint, workers=workers, threshold=threshold,
                                        min_size=1, parallel_mode='frames')
                times[name]['identification_seconds'] += time.perf_counter() - t
                t = time.perf_counter()
                tracked[name], graphs[name], states[name] = track_with_graph(labels[name], data,
                    state=states[name], return_state=True, start_time=offset,
                    sequence_id=sequence + '_' + name, **kwargs)
                times[name]['tracking_seconds'] += time.perf_counter() - t
                t = time.perf_counter()
                part = output / name / f'chunk_{offset:06d}'
                part.mkdir()
                write_graph(graphs[name], part)
                save_tracking_state(states[name], part / 'state.json')
                states[name] = load_tracking_state(part / 'state.json')
                dump(part / 'metadata.json', dict(timestamps=meta['timestamps'][offset:offset+len(data)],
                    dt_hours=1., normalized_units='mm/h', grid_attributes=grid, crop_origin_yx=[0, 0]))
                dump(part / 'frame_hashes.json', {key: [hashlib.sha256(memoryview(f)).hexdigest()
                    for f in array] for key, array in (
                        ('identified_frame_sha256', labels[name]), ('tracked_frame_sha256', tracked[name]))})
                times[name]['output_seconds'] += time.perf_counter() - t
            by_frame = defaultdict(dict)
            for node in graphs['A_1mm'].objects:
                by_frame[node.time][node.object.label] = node
            for i, frame in enumerate(data):
                absolute = offset + i
                t = time.perf_counter()
                a, b = labels['A_1mm'][i], labels['B_0p5mm'][i]
                core, partition = identify_envelope(frame, 1., .5, mode='partition', morph_structure=footprint)
                core2, system = identify_envelope(frame, 1., .5, mode='system', morph_structure=footprint)
                if not np.array_equal(core, a) or not np.array_equal(core2, a):
                    raise AssertionError('Envelope seed differs from A reference')
                row = accounting(frame, a, b, partition, system)
                links = core_memberships(a, b, absolute)
                system_links = core_memberships(a, system, absolute)
                counts = Counter(r['lower_label'] for r in links)
                sys_counts = Counter(r['lower_label'] for r in system_links)
                row.update(frame=absolute, timestamp=meta['timestamps'][absolute],
                    a_objects=int(a.max()), b_objects=int(b.max()),
                    b_objects_with_multiple_a_seeds=sum(n > 1 for n in counts.values()),
                    b_objects_without_a_seed=int(b.max()) - len(counts),
                    seeded_weak_components=len(sys_counts),
                    multi_core_weak_components=sum(n > 1 for n in sys_counts.values()))
                mappings.extend(links)
                bridge_rows.extend(dict(frame=absolute, weak_component_label=label,
                    core_labels=';'.join(str(r['core_label']) for r in system_links if r['lower_label'] == label))
                    for label, n in sys_counts.items() if n > 1)
                geometries.extend(envelope_geometry(frame, a, partition, by_frame[absolute], absolute))
                frame_rows.append(row)
                c_seconds += time.perf_counter() - t
                if plots and absolute in selected:
                    lookup = np.zeros(int(a.max()) + 1, dtype=np.int64)
                    for label, node in by_frame[absolute].items():
                        lookup[label] = node.branch_id
                    render(absolute, meta['timestamps'][absolute], frame, tracked['A_1mm'][i],
                           tracked['B_0p5mm'][i], lookup[partition], output)
                progress.update(absolute+1, phase='A/B tracked; C extent/accounting checked')
            del data, labels, tracked, graphs
        for name in modes:
            (output / name / 'SUCCESS').write_text('Experimental stream completed, not accuracy certification\n')
        audit_a, objects_a, edges_a, _, _ = inspect(output / 'A_1mm')
        audit_b, objects_b, edges_b, _, _ = inspect(output / 'B_0p5mm')
        projected, changes = projected_edges(objects_a, edges_a, objects_b, edges_b, mappings)
        equality = None
        if reference is not None:
            equality = reference_equivalence(output / 'A_1mm', reference)
        totals = {k: sum(row[k] for row in frame_rows) for k in (
            'all_provided_rain', 'a_rain', 'b_rain', 'c_rain', 'a_cells', 'b_cells', 'c_cells',
            'unassigned_lower_rain', 'unassigned_lower_cells', 'b_objects_with_multiple_a_seeds',
            'b_objects_without_a_seed', 'multi_core_weak_components')}
        for mode in ('a', 'b', 'c'):
            totals[mode + '_rain_fraction'] = totals[mode + '_rain'] / totals['all_provided_rain'] if totals['all_provided_rain'] else None
        if sha_file(source) != source_sha:
            raise ValueError('Prepared input changed during experiment')
        drift = [r['centroid_drift_cells'] for r in geometries]
        summary = dict(frames=len(rain), grid_shape=list(rain.shape[1:]),
            first_interval_end=meta['timestamps'][0], last_interval_end=meta['timestamps'][-1],
            a_structure=audit_a, b_structure=audit_b, accounting=totals,
            a_catalog_metrics=catalog_metrics(objects_a, edges_a, len(rain)),
            b_catalog_metrics=catalog_metrics(objects_b, edges_b, len(rain)),
            projected_edge_comparison=projected, a_reference_reproduction=equality,
            c_tracking_graph='Exactly reused A; no envelope tracking was run',
            c_centroid_drift_cells_quantiles={str(q): float(np.quantile(drift, q)) if drift else None
                                            for q in (.5, .9, .99, 1.)},
            stage_seconds=times, c_envelope_and_accounting_seconds=c_seconds,
            total_wall_seconds=time.perf_counter()-started, limitations=LIMITATIONS)
        write_rows(output / 'frame_accounting.csv', frame_rows)
        write_rows(output / 'a_to_b_memberships.csv', mappings)
        write_rows(output / 'c_envelope_objects.csv', geometries)
        write_rows(output / 'weak_bridge_components.csv', bridge_rows)
        dump(output / 'projected_edge_changes.json', changes)
        dump(output / 'summary.json', summary)
        (output / 'REPORT.md').write_text('# Paired 1.0 / 0.5 rain sensitivity\n\n' + LIMITATIONS +
            '\n\n```json\n' + json.dumps(summary, indent=2) + '\n```\n')
        (output / 'SUCCESS').write_text('Paired experiment and accounting completed; not accuracy validation\n')
        succeeded = True
        return summary
    finally:
        progress.finish(succeeded)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input', type=Path, help='Prepared mm/h rate NPY, not cumulative RAINNC')
    p.add_argument('--metadata', type=Path, required=True, help='Matching preparation timestamps/grid JSON')
    p.add_argument('--units', choices=['mm/h'], required=True, help='Explicit prepared-rate units confirmation')
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--preset', type=Path, default=DEFAULT_PRESET)
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--chunk-frames', type=int, default=6)
    p.add_argument('--reference', type=Path, help='Frozen A stream for exact reproduction checks')
    p.add_argument('--no-plots', action='store_true')
    a = p.parse_args()
    print(json.dumps(run(a.input, a.metadata, a.output_dir, a.preset,
        a.workers, a.chunk_frames, a.reference, not a.no_plots), indent=2))
