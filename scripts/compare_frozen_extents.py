#!/usr/bin/env python3
"""Compare 0.5/0.1-mm/h extents of an immutable saved 1-mm/h graph.

Accepts A_1mm from compare_rain_thresholds, with its parent configuration.
Replays seed identification and verifies hashes, input identity and catalog.
No tracking call, new graph, stronger-core criterion or distance cap is used.
"""
import argparse
from collections import defaultdict
import csv
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from step.envelope import expand_seeded_envelope
from step.identification import identify, identify_frame
from audit_reference_stream import inspect, input_identity
from compare_rain_thresholds import sha_file, dump, write_rows, core_memberships
from frame_progress import FrameProgress
from run_npy_validation import disk

ROOT = Path(__file__).resolve().parents[1]
THRESHOLDS = (.5, .1)
NOTE = ('Extent sensitivity only, not storm accuracy or convection classification. '
        'Exactly one saved graph is referenced; it is not retracked. Coreless weak rain '
        'is unassigned and seed-free lifecycle stages are not recovered. Multiple seeds '
        'in a weak component are a grouping diagnostic, not a false-merge verdict. '
        'Rain sums use equal grid cells and supplied RAINNC only. Distances are grid '
        'cells to ANY seed, not owned-seed or wet-path distances. No extent cap is applied. '
        'Nonseed ownership is checked separately from seed identity and support.')


def digest_frame(labels):
    return hashlib.sha256(memoryview(np.ascontiguousarray(labels))).hexdigest()


def manifest(reference):
    files = [reference.parent / 'configuration.json', reference / 'SUCCESS']
    files.extend(reference.parent / name
                 for name in ('frame_accounting.csv', 'c_envelope_objects.csv')
                 if (reference.parent / name).exists())
    files.extend(p for part in sorted(reference.glob('chunk_*'))
                 for p in sorted(part.iterdir()) if p.is_file())
    return {str(p.resolve()): sha_file(p) for p in files}


def extent_frame(rain, cores, nodes, frame, threshold, distance):
    valid = np.isfinite(rain)
    wet = valid & (rain >= threshold)
    components, count = ndimage.label(wet, np.ones((3, 3), dtype=bool))
    membership = core_memberships(cores, components, frame)
    seeded_ids = sorted({r['lower_label'] for r in membership})
    seeded = np.zeros(count + 1, dtype=bool)
    seeded[seeded_ids] = True
    expected = seeded[components]
    partition = expand_seeded_envelope(rain, cores, threshold)
    if not np.array_equal(partition > 0, expected):
        raise AssertionError('Watershed support differs from seed-containing wet components')
    if not np.array_equal(partition[cores > 0], cores[cores > 0]):
        raise AssertionError('Seed identity was changed')
    if np.any((partition > 0) & ~wet):
        raise AssertionError('Extent crossed dry/missing pixels')
    if set(np.unique(partition[partition > 0])) != set(nodes):
        raise AssertionError('Extent IDs differ from the frame catalog')
    n = int(cores.max()) + 1
    active = np.flatnonzero(partition)
    labels = partition.ravel()[active]
    weights = rain.ravel()[active].astype(np.float64)
    size = np.bincount(labels, minlength=n)
    mass = np.bincount(labels, weights=weights, minlength=n)
    sy = np.bincount(labels, weights=weights * (active // rain.shape[1]), minlength=n)
    sx = np.bincount(labels, weights=weights * (active % rain.shape[1]), minlength=n)
    radii = ndimage.maximum(distance, labels=partition, index=list(nodes)) if nodes else []
    geometry = []
    for label, radius in zip(nodes, radii):
        node = nodes[label]
        y, x = sy[label] / mass[label], sx[label] / mass[label]
        geometry.append(dict(frame=frame, threshold_mm_h=threshold, core_label=label,
            core_node_id=int(node['node_id']), branch_id=int(node['branch_id']),
            core_cells=int(node['area_cells']), envelope_cells=int(size[label]),
            envelope_rain=float(mass[label]),
            envelope_to_core_cell_ratio=float(size[label] / int(node['area_cells'])),
            centroid_drift_cells=float(np.hypot(y-float(node['centroid_y']), x-float(node['centroid_x']))),
            farthest_distance_to_any_seed_cells=float(radius)))
    total = lambda mask: float(np.sum(rain[mask], dtype=np.float64))
    assigned, unassigned = partition > 0, wet & (partition == 0)
    stats = dict(frame=frame, threshold_mm_h=threshold,
        all_provided_rain=total(valid), core_cells=int((cores > 0).sum()),
        core_rain=total(cores > 0), eligible_cells=int(wet.sum()),
        eligible_rain=total(wet), assigned_cells=int(assigned.sum()),
        assigned_rain=total(assigned), unassigned_cells=int(unassigned.sum()),
        unassigned_rain=total(unassigned), wet_components=count,
        seeded_weak_components=len(seeded_ids), unseeded_weak_components=count-len(seeded_ids))
    if not np.isclose(stats['assigned_rain'] + stats['unassigned_rain'],
                      stats['eligible_rain'], rtol=1e-12, atol=1e-8):
        raise AssertionError('Assigned/unassigned eligible rain is not conserved')
    groups = defaultdict(list)
    for row in membership:
        groups[row['lower_label']].append(row['core_label'])
    stats['multi_seed_weak_components'] = sum(len(ids) > 1 for ids in groups.values())
    bridges = [dict(frame=frame, threshold_mm_h=threshold,
        weak_component_label=component, seed_count=len(ids),
        core_labels=';'.join(map(str, sorted(ids))))
        for component, ids in groups.items() if len(ids) > 1]
    return partition, stats, geometry, bridges


def support_change(rain, cores, high, low, frame):
    if np.any((high > 0) & (low == 0)):
        raise AssertionError('Lower threshold lost assigned support')
    if not np.array_equal(low[cores > 0], cores[cores > 0]):
        raise AssertionError('Lower threshold changed seed IDs')
    added = (low > 0) & (high == 0)
    changed = (high > 0) & (low != high)
    # Added support may include old >=0.5 rain newly connected by a weak bridge.
    return dict(frame=frame, added_cells=int(added.sum()),
        added_rain=float(np.sum(rain[added], dtype=np.float64)),
        added_below_0p5_cells=int((added & (rain < .5)).sum()),
        added_existing_ge_0p5_cells=int((added & (rain >= .5)).sum()),
        added_existing_ge_0p5_rain=float(np.sum(rain[added & (rain >= .5)], dtype=np.float64)),
        reassigned_existing_cells=int(changed.sum()),
        reassigned_existing_rain=float(np.sum(rain[changed], dtype=np.float64)))


def quantiles(values):
    return {str(q): float(np.quantile(values, q)) if values else None
            for q in (.5, .9, .99, 1.)}


def previous_extent_equivalence(reference, frame_rows, geometry):
    """Compare saved C accounting/geometry, not unavailable old C raster hashes."""
    accounting_path = reference.parent / 'frame_accounting.csv'
    geometry_path = reference.parent / 'c_envelope_objects.csv'
    if not accounting_path.exists() or not geometry_path.exists():
        return None
    with accounting_path.open(newline='') as handle:
        old_frames = list(csv.DictReader(handle))
    with geometry_path.open(newline='') as handle:
        old_objects = list(csv.DictReader(handle))
    new_frames = [r for r in frame_rows if r['threshold_mm_h'] == .5]
    new_objects = [r for r in geometry if r['threshold_mm_h'] == .5]
    if len(old_frames) != len(new_frames) or len(old_objects) != len(new_objects):
        raise AssertionError('Previous 0.5 extent row counts differ')
    for old, new in zip(old_frames, new_frames):
        for old_key, new_key in (('frame', 'frame'), ('c_cells', 'assigned_cells'),
                                 ('c_rain', 'assigned_rain'), ('a_rain', 'core_rain'),
                                 ('unassigned_lower_cells', 'unassigned_cells'),
                                 ('unassigned_lower_rain', 'unassigned_rain')):
            if not np.isclose(float(old[old_key]), new[new_key], rtol=1e-12, atol=1e-8):
                raise AssertionError(f'Previous 0.5 accounting differs: {old_key}')
    # Catalog rows can have different ordering; compare immutable frame/label keys.
    old_by_key = {(int(r['frame']), int(r['core_label'])): r for r in old_objects}
    if len(old_by_key) != len(old_objects):
        raise AssertionError('Previous 0.5 geometry has duplicate identities')
    for new in new_objects:
        old = old_by_key[(new['frame'], new['core_label'])]
        for key in old:
            if not np.isclose(float(old[key]), new[key], rtol=1e-12, atol=1e-8):
                raise AssertionError(f'Previous 0.5 geometry differs: {key}')
    return dict(accounting_equal=True, geometry_equal=True,
                note='Old C raster hashes were not saved; this is not pixelwise ownership equivalence.')


def render(frame, stamp, rain, cores, nodes, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), layout='constrained')
    im = axes[0].imshow(rain, origin='lower', cmap='Blues', vmin=0, vmax=10)
    fig.colorbar(im, ax=axes[0], label='mm/h')
    axes[0].set_title('Rain; display capped at 10')
    lookup = np.zeros(int(cores.max())+1, dtype=np.int64)
    for label, node in nodes.items():
        lookup[label] = int(node['branch_id'])
    for ax, threshold in zip(axes[1:], THRESHOLDS):
        labels = expand_seeded_envelope(rain, cores, threshold)
        branches = lookup[labels]
        ax.imshow(np.ma.masked_where(branches == 0, branches % 20), origin='lower',
                  cmap='tab20', vmin=0, vmax=19, interpolation='nearest')
        ax.contour(cores > 0, levels=[.5], colors='black', linewidths=.25)
        ax.set_title(f'{threshold:g} extent; same frozen branch IDs')
    for ax in axes:
        ax.set(xlabel='x (grid cells)', ylabel='y (grid cells)')
    fig.suptitle(f'{stamp} | frame {frame}; black: 1-mm/h seed support; colors repeat')
    fig.savefig(output / f'extent_{frame:03d}.png', dpi=140)
    plt.close(fig)


def run(source, metadata_path, reference, output, workers=1, chunk_frames=6, plots=True):
    if output.exists():
        raise FileExistsError(output)
    if workers < 1 or chunk_frames < 1:
        raise ValueError('workers and chunk-frames must be positive')
    started = time.perf_counter()
    config_path = reference.parent / 'configuration.json'
    config = json.loads(config_path.read_text())
    preset = config['preset']
    if (preset['identification']['threshold_mm_h'] != 1 or
            preset['input_contract'] != dict(grid_km=4, dt_hours=1,
                rain_units_after_normalization='mm/h')):
        raise ValueError('Requires an explicitly prepared 4-km hourly 1-mm/h reference')
    source_sha, metadata_sha = sha_file(source), sha_file(metadata_path)
    if (source_sha != config['input_sha256'] or metadata_sha != config['metadata_sha256']):
        raise ValueError('Prepared input / metadata checksum differs from reference')
    rain = np.load(source, mmap_mode='r')
    meta = json.loads(metadata_path.read_text())
    stamps = np.asarray(meta['timestamps'], dtype='datetime64[s]')
    if (rain.ndim != 3 or not len(rain) or not np.issubdtype(rain.dtype, np.floating) or
            len(stamps) != len(rain) or np.isnat(stamps).any() or
            not np.all(np.diff(stamps) == np.timedelta64(1, 'h'))):
        raise ValueError('Requires floating hourly rate cube and complete hourly chronology')
    structure, objects, _, _, _ = inspect(reference)
    identity = input_identity(reference)
    attributes = dict(identity['grid']['grid_attributes'])
    attributes.pop('shape', None)
    expected_attributes = dict(meta['grid'])
    expected_attributes.pop('shape', None)
    if (identity['timestamps'] != meta['timestamps'] or structure['frames'] != len(rain) or
            list(rain.shape[1:]) != meta['grid']['shape'] or
            identity['grid']['processed_grid_shape'] != meta['grid']['shape'] or
            attributes != expected_attributes or meta['grid']['DX'] != 4000 or
            meta['grid']['DY'] != 4000 or identity['grid']['dt_hours'] != 1 or
            identity['grid']['normalized_units'] != 'mm/h' or
            identity['grid']['crop_origin_yx'] != [0, 0]):
        raise ValueError('Saved reference chronology / grid / units differ')
    reference_hashes = [h for part in sorted(reference.glob('chunk_*'))
        for h in json.loads((part / 'frame_hashes.json').read_text())['identified_frame_sha256']]
    if len(reference_hashes) != len(rain):
        raise ValueError('Reference does not provide every identified-frame hash')
    before = manifest(reference)
    nodes = defaultdict(dict)
    for row in objects:
        nodes[int(row['time'])][int(row['local_label'])] = row
    output.mkdir(parents=True)
    dump(output / 'configuration.json', dict(prepared_input=str(source.resolve()),
        metadata=str(metadata_path.resolve()), reference=str(reference.resolve()),
        input_sha256=source_sha, metadata_sha256=metadata_sha, preset=preset,
        thresholds_mm_h=list(THRESHOLDS), connectivity=8, distance_cap=None,
        workers=workers, chunk_frames=chunk_frames, source_graph_files_sha256=before,
        python_version=sys.version, dependency_versions={name: version(name)
            for name in ('numpy', 'scipy', 'scikit-image')},
        code_sha256={name: sha_file(ROOT / name) for name in
            ('scripts/compare_frozen_extents.py', 'step/envelope.py', 'step/identification.py',
             'step/tracking.py')},
        git_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        limitations=NOTE))
    progress = FrameProgress(len(rain))
    progress.unit, progress.verb = 'hours', 'audited'
    frame_rows, geometry, bridges, differences, seed_hashes = [], [], [], [], []
    id_seconds, extent_seconds = 0., 0.
    footprint = disk(preset['identification']['bridge_radius_cells'])
    succeeded = False
    try:
        for offset in range(0, len(rain), chunk_frames):
            progress.update(offset, phase='replaying frozen seeds; no tracking changes')
            data = np.array(rain[offset:offset+chunk_frames], copy=True)
            if np.any(data[np.isfinite(data)] < 0):
                raise ValueError('Negative hourly rate; cumulative resets must be explicit')
            t = time.perf_counter()
            core_cube = identify(data, footprint, workers=workers,
                threshold=1, min_size=preset['identification']['min_size_cells'])
            id_seconds += time.perf_counter()-t
            for i, current in enumerate(data):
                frame, cores = offset+i, core_cube[i]
                seed_hash = digest_frame(cores)
                if seed_hash != reference_hashes[frame]:
                    raise AssertionError(f'Frozen seed raster differs at frame {frame}')
                if set(np.unique(cores[cores > 0])) != set(nodes[frame]):
                    raise AssertionError('Seed labels differ from frozen object catalog')
                counts = np.bincount(cores.ravel())
                if any(counts[label] != int(node['area_cells'])
                       for label, node in nodes[frame].items()):
                    raise AssertionError('Seed areas differ from frozen catalog')
                seed_hashes.append(seed_hash)
                t = time.perf_counter()
                distance = ndimage.distance_transform_edt(cores == 0) if cores.any() else np.zeros(cores.shape)
                partitions = []
                for threshold in THRESHOLDS:
                    labels, stats, geom, groups = extent_frame(current, cores, nodes[frame],
                                                             frame, threshold, distance)
                    stats['timestamp'] = meta['timestamps'][frame]
                    frame_rows.append(stats)
                    geometry.extend(geom)
                    bridges.extend(groups)
                    partitions.append(labels)
                differences.append(support_change(current, cores, *partitions, frame))
                extent_seconds += time.perf_counter()-t
                progress.update(frame+1, phase='0.5/0.1 extents checked; frozen graph reused')
            del data, core_cube, partitions, distance, labels
        totals = {}
        numeric = ('all_provided_rain', 'core_rain', 'core_cells', 'eligible_rain',
                   'eligible_cells', 'assigned_rain', 'assigned_cells',
                   'unassigned_rain', 'unassigned_cells', 'wet_components',
                   'seeded_weak_components', 'unseeded_weak_components',
                   'multi_seed_weak_components')
        for threshold in THRESHOLDS:
            rows = [r for r in frame_rows if r['threshold_mm_h'] == threshold]
            values = {key: sum(r[key] for r in rows) for key in numeric}
            denominator = values['all_provided_rain']
            values['assigned_rain_fraction'] = values['assigned_rain']/denominator if denominator else None
            values['eligible_rain_fraction'] = values['eligible_rain']/denominator if denominator else None
            geom = [r for r in geometry if r['threshold_mm_h'] == threshold]
            values['centroid_drift_cells_quantiles'] = quantiles([r['centroid_drift_cells'] for r in geom])
            values['area_ratio_quantiles'] = quantiles([r['envelope_to_core_cell_ratio'] for r in geom])
            values['max_any_seed_distance_cells_quantiles'] = quantiles(
                [r['farthest_distance_to_any_seed_cells'] for r in geom])
            values['largest_seed_count_in_weak_component'] = max(
                (r['seed_count'] for r in bridges if r['threshold_mm_h'] == threshold), default=1 if geom else 0)
            totals[str(threshold)] = values
        change = {key: sum(r[key] for r in differences) for key in differences[0] if key != 'frame'}
        extremes = {name: sorted(geometry, key=lambda r: r[key], reverse=True)[:10]
            for name, key in (('area_ratio', 'envelope_to_core_cell_ratio'),
                              ('centroid_drift', 'centroid_drift_cells'),
                              ('any_seed_distance', 'farthest_distance_to_any_seed_cells'))}
        if plots:
            progress.update(phase='finalizing diagnostic maps; not yet successful')
            selected = set(np.linspace(0, len(rain)-1, min(4, len(rain)), dtype=int).tolist())
            selected.update(rows[0]['frame'] for rows in extremes.values() if rows)
            for frame in sorted(selected):
                current = np.array(rain[frame], copy=True)
                cores = identify_frame(current, footprint,
                    min_size=preset['identification']['min_size_cells'], threshold=1)
                if digest_frame(cores) != reference_hashes[frame]:
                    raise AssertionError('Plot seed replay differs')
                render(frame, meta['timestamps'][frame], current, cores, nodes[frame], output)
        previous = previous_extent_equivalence(reference, frame_rows, geometry)
        if (sha_file(source) != source_sha or sha_file(metadata_path) != metadata_sha or
                manifest(reference) != before):
            raise AssertionError('Input, metadata or saved reference changed during experiment')
        summary = dict(frames=len(rain), grid_shape=list(rain.shape[1:]),
            first_interval_end=meta['timestamps'][0], last_interval_end=meta['timestamps'][-1],
            tracking_graph='Reused immutable saved A graph; not rerun', reference_structure=structure,
            exact_seed_hashes_equal=True, graph_source_files_unchanged=True,
            catalog_seed_ids_and_areas_equal=True, accounting_and_barrier_checks_pass=True,
            previous_0p5_extent_comparison=previous,
            extents=totals, lower_threshold_support_changes=change,
            identification_seconds=id_seconds, extent_and_statistics_seconds=extent_seconds,
            total_wall_seconds=time.perf_counter()-started, limitations=NOTE)
        write_rows(output / 'frame_accounting.csv', frame_rows)
        write_rows(output / 'extent_objects.csv', geometry)
        write_rows(output / 'weak_bridge_components.csv', bridges)
        write_rows(output / 'support_changes.csv', differences)
        dump(output / 'seed_frame_hashes.json', seed_hashes)
        dump(output / 'extreme_objects.json', extremes)
        dump(output / 'summary.json', summary)
        (output / 'REPORT.md').write_text('# Frozen tracking: 0.5 / 0.1 rain extents\n\n' + NOTE +
            '\n\n```json\n' + json.dumps(summary, indent=2) + '\n```\n')
        (output / 'SUCCESS').write_text('Extent sensitivity checks completed; no threshold promotion\n')
        succeeded = True
        return summary
    finally:
        progress.finish(succeeded)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, help='Prepared hourly mm/h NPY, not cumulative RAINNC')
    parser.add_argument('--metadata', required=True, type=Path)
    parser.add_argument('--reference', required=True, type=Path, help='Paired A_1mm output directory')
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--units', required=True, choices=['mm/h'])
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--chunk-frames', type=int, default=6)
    parser.add_argument('--no-plots', action='store_true')
    args = parser.parse_args()
    print(json.dumps(run(args.input, args.metadata, args.reference, args.output_dir,
                         args.workers, args.chunk_frames, not args.no_plots), indent=2))
