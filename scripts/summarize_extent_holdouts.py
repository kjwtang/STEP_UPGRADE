"""Describe successful fixed-seed extent runs without selecting a new threshold."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from compare_rain_thresholds import sha_file

CHECKS = ('exact_seed_hashes_equal', 'graph_source_files_unchanged',
          'catalog_seed_ids_and_areas_equal', 'accounting_and_barrier_checks_pass')
BINS = (('1 cell', 1, 1), ('2-15 cells', 2, 15),
        ('16-999 cells', 16, 999), ('>=1000 cells', 1000, float('inf')))


def read_run(path):
    if not (path / 'SUCCESS').is_file():
        raise ValueError(f'Incomplete run: {path}')
    summary = json.loads((path / 'summary.json').read_text())
    if not all(summary.get(key) is True for key in CHECKS):
        raise ValueError(f'Integrity checks did not pass: {path}')
    config = json.loads((path / 'configuration.json').read_text())
    with (path / 'extent_objects.csv').open(newline='') as handle:
        objects = list(csv.DictReader(handle))
    groups = []
    for threshold, total in summary['extents'].items():
        rows = [r for r in objects if float(r['threshold_mm_h']) == float(threshold)]
        rain = sum(float(r['envelope_rain']) for r in rows)
        if not np.isclose(rain, total['assigned_rain'], rtol=1e-12, atol=1e-8):
            raise AssertionError('Object rain does not sum to assigned rain')
        for label, low, high in BINS:
            selected = [r for r in rows if low <= int(r['core_cells']) <= high]
            def quantile(key, q):
                return float(np.quantile([float(r[key]) for r in selected], q)) if selected else None
            groups.append(dict(threshold_mm_h=float(threshold), seed_area_bin=label,
                object_snapshots=len(selected),
                assigned_rain_fraction=sum(float(r['envelope_rain']) for r in selected)/rain if rain else None,
                area_ratio_p99=quantile('envelope_to_core_cell_ratio', .99),
                area_ratio_max=quantile('envelope_to_core_cell_ratio', 1),
                centroid_drift_p99_cells=quantile('centroid_drift_cells', .99)))
    return dict(path=str(path.resolve()), summary=summary, seed_area_bins=groups,
        preset=config['preset'], artifacts_sha256={name: sha_file(path / name) for name in
            ('SUCCESS', 'summary.json', 'configuration.json', 'extent_objects.csv')})


def run(paths, output):
    if output.exists():
        raise FileExistsError(output)
    if not paths:
        raise ValueError('At least one successful extent run is required')
    runs = [read_run(path) for path in paths]
    if any(r['summary']['extents']['0.1']['all_provided_rain'] <= 0 for r in runs):
        raise ValueError('Rain-fraction report requires a positive provided-rain total')
    if any(r['preset'] != runs[0]['preset'] or
           r['summary']['grid_shape'] != runs[0]['summary']['grid_shape'] or
           r['summary']['frames'] != runs[0]['summary']['frames'] for r in runs):
        raise ValueError('Compare matched preset, grid and duration only')
    starts = [r['summary']['first_interval_end'] for r in runs]
    if len(set(starts)) != len(starts):
        raise ValueError('Duplicate date windows are not cross-date validation')
    lines = ['# Fixed-seed extent cross-date validation', '',
        'Same committed preset, grid and duration. Seed identity, graph source integrity,',
        'catalog areas, barriers and rain conservation passed in every supplied run.',
        'Different dates in one model/year are not independent observations or tracking truth.', '',
        '| First interval end | Extent mm/h | Assigned rain % | Unassigned eligible rain % of all rain | Multi-seed components | Max seed labels/component | Area ratio p99 / max | Drift p99 cells |',
        '|---|---:|---:|---:|---:|---:|---:|---:|']
    for record in runs:
        summary = record['summary']
        for threshold, stats in summary['extents'].items():
            ratio = stats['area_ratio_quantiles']
            lines.append(f"| {summary['first_interval_end']} | {threshold} | "
                f"{100*stats['assigned_rain_fraction']:.3f} | "
                f"{100*stats['unassigned_rain']/stats['all_provided_rain']:.3f} | "
                f"{stats['multi_seed_weak_components']} | {stats['largest_seed_count_in_weak_component']} | "
                f"{ratio['0.99']:.2f} / {ratio['1.0']:.0f} | "
                f"{stats['centroid_drift_cells_quantiles']['0.99']:.2f} |")
    lines += ['', '## Lower-threshold support', '',
        '| First interval end | Extra rain percentage points | Previously coreless >=0.5 pixels added | Their % of extra rain | Existing assigned pixel ownership changes |',
        '|---|---:|---:|---:|---:|']
    for record in runs:
        s = record['summary']
        delta = s['lower_threshold_support_changes']
        extra = delta['added_rain']
        lines.append(f"| {s['first_interval_end']} | "
            f"{100*(s['extents']['0.1']['assigned_rain_fraction']-s['extents']['0.5']['assigned_rain_fraction']):.3f} | "
            f"{delta['added_existing_ge_0p5_cells']} | "
            f"{100*delta['added_existing_ge_0p5_rain']/extra if extra else 0:.2f} | "
            f"{delta['reassigned_existing_cells']} |")
    lines += ['', '## Descriptive seed-area subsets at 0.1 mm/h', '',
        'Bins describe existing seed snapshots; they are not a new core-area filter.', '',
        '| First interval end | Seed cells | Snapshots | % assigned rain | Area ratio p99 / max | Drift p99 cells |',
        '|---|---|---:|---:|---:|---:|']
    for record in runs:
        for group in record['seed_area_bins']:
            if group['threshold_mm_h'] != .1 or not group['object_snapshots']:
                continue
            lines.append(f"| {record['summary']['first_interval_end']} | {group['seed_area_bin']} | "
                f"{group['object_snapshots']} | {100*group['assigned_rain_fraction']:.3f} | "
                f"{group['area_ratio_p99']:.2f} / {group['area_ratio_max']:.2f} | "
                f"{group['centroid_drift_p99_cells']:.2f} |")
    lines += ['', '## Interpretation limits', '',
        'Coverage is not accuracy. Multi-seed components do not certify a false merge.',
        'Centroid drift is core-to-extent measurement drift, not motion. Rain sums use',
        'equal grid weights and provided RAINNC only; weak coreless lifecycle phases',
        'remain excluded. Expanded centroids never enter the frozen tracker.',
        'No default, seed requirement, distance cap or tracking rule is promoted here.', '',
        'Per-object bins and source artifact hashes are recorded in `summary.json`.', '']
    output.mkdir(parents=True)
    (output / 'summary.json').write_text(json.dumps(dict(runs=runs), indent=2)+'\n')
    (output / 'REPORT.md').write_text('\n'.join(lines))
    return runs


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs', nargs='+', type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    run(args.runs, args.output_dir)
    print(f'Report: {args.output_dir / "REPORT.md"}')
