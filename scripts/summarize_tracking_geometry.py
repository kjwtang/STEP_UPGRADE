"""Compare experimental graphs using identical authoritative seed-node weights."""
import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np

from audit_reference_stream import inspect, input_identity
from compare_frozen_extents import manifest
from compare_rain_thresholds import dump, sha_file, reference_equivalence
from compare_tracking_geometry import GEOMETRIES, NOTE, compare_links

AREA_BINS = ((1, 1, '1'), (2, 15, '2-15'), (16, 999, '16-999'),
             (1000, float('inf'), '>=1000'))


def event_seed_sizes(seed_objects, events):
    seed = {int(r['node_id']): r for r in seed_objects}
    counts = {name: 0 for _, _, name in AREA_BINS}
    for event in events:
        ids = [int(n) for k in ('parent_node_ids', 'child_node_ids')
               for n in event[k].split(';') if n]
        minimum = min(int(seed[n]['area_cells']) for n in ids)
        for low, high, name in AREA_BINS:
            counts[name] += low <= minimum <= high
    return dict(events_by_smallest_original_seed_cells=counts,
                note='Descriptive endpoint sizes; not a new event filter or accuracy label.')


def edge_change_seed_sizes(seed_objects, changes):
    seed = {(int(r['time']), int(r['local_label'])): r for r in seed_objects}
    bins = {name: 0 for _, _, name in AREA_BINS}
    for change in changes:
        p, c = tuple(change['pair'][:2]), tuple(change['pair'][2:])
        minimum = min(int(seed[k]['area_cells']) for k in (p, c))
        for low, high, name in AREA_BINS:
            bins[name] += low <= minimum <= high
    return bins


def fixed_seed_metrics(seed_objects, objects, edges, frames):
    seed = {int(r['node_id']): r for r in seed_objects}
    current = {int(r['node_id']): r for r in objects}
    if {n: (r['time'], r['local_label']) for n, r in seed.items()} != {
            n: (r['time'], r['local_label']) for n, r in current.items()}:
        raise ValueError('Seed identity changed')
    incoming = {int(r['child_node_id']) for r in edges}
    eligible = [n for n, r in seed.items() if int(r['time']) > 0]
    associated = [n for n in eligible if n in incoming]
    mass = lambda ids: sum(float(seed[n]['precipitation_sum']) for n in ids)
    denominator = mass(eligible)
    groups = defaultdict(list)
    for r in objects:
        groups[int(r['branch_id'])].append(int(r['time']))
    spans = [max(v)-min(v)+1 for v in groups.values()]
    return dict(incoming_seed_node_fraction=len(associated)/len(eligible) if eligible else None,
        incoming_seed_rain_fraction=mass(associated)/denominator if denominator else None,
        seed_rain_denominator=denominator,
        one_snapshot_branch_fraction=sum(len(v) == 1 for v in groups.values())/len(groups) if groups else None,
        branch_span_hours_median=float(np.median(spans)) if spans else None,
        branch_span_hours_max=max(spans, default=None),
        branches_touching_window=sum(0 in v or frames-1 in v for v in groups.values()),
        note='All rain weights use the same original >=1-mm/h catalog. Association is not accuracy; '
             'branch spans include gaps, event cuts and window censoring.')


def run(paths, output):
    if output.exists():
        raise FileExistsError(output)
    records = []
    for path in paths:
        if not (path / 'SUCCESS').is_file():
            raise ValueError(f'Incomplete experiment: {path}')
        summary = json.loads((path / 'summary.json').read_text())
        if not all(summary['baseline_control'].values()) or not summary['source_files_unchanged']:
            raise ValueError('Control or integrity checks failed')
        plan = json.loads((path / 'PLAN.json').read_text())
        extents = Path(plan['source'])
        config = json.loads((extents / 'configuration.json').read_text())
        reference = Path(config['reference'])
        if manifest(reference) != plan['source_graph_files_sha256']:
            raise ValueError('Original graph changed after experiment')
        if {n: sha_file(extents / n) for n in plan['extent_source_sha256']} != plan['extent_source_sha256']:
            raise ValueError('Extent source changed after experiment')
        if sha_file(Path(config['prepared_input'])) != plan['input_sha256']:
            raise ValueError('Prepared input changed after experiment')
        if sha_file(Path(config['metadata'])) != plan['metadata_sha256']:
            raise ValueError('Prepared metadata changed after experiment')
        reference_equivalence(path / 'seed_control', reference)
        _, seed, _, _, _ = inspect(reference)
        variants = {}
        for name in GEOMETRIES:
            structure, objects, edges, events, _ = inspect(path / name)
            variants[name] = dict(structure=structure,
                fixed_seed_metrics=fixed_seed_metrics(seed, objects, edges, summary['frames']),
                event_seed_sizes=event_seed_sizes(seed, events))
            if name != 'seed_control':
                changes = json.loads((path / (name+'_edge_changes.json')).read_text())
                fresh, fresh_changes = compare_links(reference, path / name)
                if fresh != summary['comparisons'][name] or fresh_changes != changes:
                    raise ValueError('Saved comparison differs from current experimental graph')
                variants[name]['edge_changes_by_smallest_seed_cells'] = {
                    kind: edge_change_seed_sizes(seed, rows) for kind, rows in changes.items()}
        records.append(dict(path=str(path.resolve()), plan=plan, summary=summary, variants=variants,
            input_contract=input_identity(path / 'seed_control')['grid'],
            artifacts_sha256={n: sha_file(path / n) for n in ('PLAN.json', 'summary.json', 'SUCCESS')}))
    if not records or len({r['summary']['first_interval_end'] for r in records}) != len(records):
        raise ValueError('Require nonempty distinct date windows')
    if any(r['plan']['preset'] != records[0]['plan']['preset'] for r in records):
        raise ValueError('Tracking parameters differ')
    if any(r['summary']['frames'] != records[0]['summary']['frames'] or
           r['input_contract'] != records[0]['input_contract'] for r in records):
        raise ValueError('Require matched grid, cadence, units and duration')
    lines = ['# Fixed-seed tracking geometry sensitivity', '', NOTE, '',
        'All controls reproduced the saved complete graphs. All experiments retained every',
        'original node/frame/local-label identity and passed graph structural checks.',
        'Only experimental features and resulting lineage changed; no reference was replaced.', '',
        '| First interval | Geometry | Nodes | Edges | Gap edges | Events split / merge / complex | Branches | Associated seed rain % | One-snapshot branches % | Max branch span h |',
        '|---|---|---:|---:|---:|---|---:|---:|---:|---:|']
    for r in records:
        for name, v in r['variants'].items():
            s, m = v['structure'], v['fixed_seed_metrics']
            pct = lambda x: f'{100*x:.3f}' if x is not None else 'n/a'
            events = ' / '.join(str(s['events'].get(k, 0)) for k in ('split', 'merge', 'complex'))
            lines.append(f"| {r['summary']['first_interval_end']} | {name} | {s['nodes']} | {s['edges']} | "
                f"{s['edge_types'].get('gap_continue', 0)} | {events} | {s['branches']} | "
                f"{pct(m['incoming_seed_rain_fraction'])} | {pct(m['one_snapshot_branch_fraction'])} | "
                f"{m['branch_span_hours_max']} |")
    lines += ['', '## Edge differences relative to immutable seed tracking', '',
        '| First interval | Geometry | Retained endpoint pairs | Added | Removed | Event type changed on retained pairs |',
        '|---|---|---:|---:|---:|---:|']
    for r in records:
        for name, c in r['summary']['comparisons'].items():
            lines.append(f"| {r['summary']['first_interval_end']} | {name} | {c['retained_links']} | "
                         f"{c['added']} | {c['removed']} | {c['changed_event']} |")
    lines += ['', '## Interpretation', '',
        'Fixed seed rain denominators avoid counting newly assigned weak rain as an apparent',
        'tracking improvement. Higher association and fewer branches are still not ground truth.',
        'Raw branch numbers/colors are not cross-run identities. Changes are keyed by immutable',
        '(parent frame, local seed label, child frame, local seed label).',
        'The three feature geometries jointly change overlap, centroid, area and mean intensity.',
        'This is not a single-feature causal ablation or a calibrated threshold choice.',
        'Different dates from one model/year are not independent observational validation.',
        'Do not promote an expanded geometry without false-link/event review and independent',
        'time/year testing. Optional strong-phase sidecar labels do not enter any tracker here.', '']
    output.mkdir(parents=True)
    dump(output / 'summary.json', dict(runs=records))
    (output / 'REPORT.md').write_text('\n'.join(lines))
    dump(output / 'SUCCESS', dict(scope='Descriptive sensitivity, not tracking accuracy'))
    return records


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs', nargs='+', type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    run(args.runs, args.output_dir)
    print(f'Report: {args.output_dir / "REPORT.md"}')
