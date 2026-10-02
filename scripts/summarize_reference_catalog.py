#!/usr/bin/env python3
"""Rain-catalog coverage/size diagnostics, not storm truth or a tracking accuracy score."""
import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np
from audit_reference_stream import inspect


def summarize(source, output):
    if output.exists():
        raise FileExistsError(output)
    structural, objects, edges, events, _ = inspect(source)
    contract_path = source / 'execution_contract.json'
    settings_path = source / 'settings.json'
    if contract_path.exists():
        threshold = json.loads(contract_path.read_text())['preset']['identification']['threshold_mm_h']
    elif settings_path.exists():
        threshold = json.loads(settings_path.read_text()).get('threshold')
    else:
        threshold = None
    branches = defaultdict(list)
    families = defaultdict(list)
    areas, mass, times = [], [], []
    for row in objects:
        branches[int(row['branch_id'])].append(row)
        families[int(row['family_id'])].append(row)
        areas.append(int(row['area_cells']))
        mass.append(float(row['precipitation_sum']))
        times.append(int(row['time']))
    areas, mass, times = np.array(areas), np.array(mass), np.array(times)
    total = float(mass.sum())
    grouped = []
    for low, high in ((1, 16), (16, 100), (100, 1000), (1000, None)):
        mask = areas >= low
        if high is not None:
            mask &= areas < high
        grouped.append(dict(min_cells=low, max_cells_exclusive=high,
            objects=int(mask.sum()), object_fraction=float(mask.mean()) if len(mask) else 0.,
            labeled_rain_mm_cells=float(mass[mask].sum()),
            labeled_rain_fraction=float(mass[mask].sum() / total) if total else None))
    incoming = {int(r['child_node_id']) for r in edges}
    outgoing = {int(r['parent_node_id']) for r in edges}
    association = []
    for floor in (1, 16, 100, 1000):
        selected = [r for r in objects if int(r['area_cells']) >= floor]
        for direction, ids, eligible in (
                ('incoming', incoming, [r for r in selected if int(r['time']) > 0]),
                ('outgoing', outgoing, [r for r in selected if int(r['time']) < structural['frames'] - 1])):
            linked = [r for r in eligible if int(r['node_id']) in ids]
            denominator = sum(float(r['precipitation_sum']) for r in eligible)
            association.append(dict(min_area_cells=floor, direction=direction,
                eligible_objects=len(eligible), linked_objects=len(linked),
                object_fraction=len(linked) / len(eligible) if eligible else None,
                labeled_rain_fraction=sum(float(r['precipitation_sum']) for r in linked) / denominator
                    if denominator else None))
    spans = [max(int(r['time']) for r in rows) - min(int(r['time']) for r in rows) + 1
             for rows in branches.values()]
    largest_families = sorted((dict(family_id=family, nodes=len(rows),
        first_frame=min(int(r['time']) for r in rows), last_frame=max(int(r['time']) for r in rows),
        branches=len({r['branch_id'] for r in rows}),
        labeled_rain_mm_cells=sum(float(r['precipitation_sum']) for r in rows))
        for family, rows in families.items()), key=lambda r: r['nodes'], reverse=True)[:10]
    report = dict(structure=structural, identification_threshold_mm_h=threshold,
        area_groups=grouped, association_coverage=association,
        branch_span_hours_median=float(np.median(spans)) if spans else None,
        branch_span_hours_max=max(spans, default=None),
        one_snapshot_branch_fraction=sum(len(v) == 1 for v in branches.values()) / len(branches)
            if branches else None,
        largest_final_families=largest_families,
        limitations='Only labeled rain under the saved identification contract, not total rainfall. mm-cells sums use equal grid-cell weighting, not geographic area-corrected rain volume. Associations are accepted graph links, not accuracy. Domain/time boundaries censor tracks; no observational or convective classification is implied.')
    output.mkdir(parents=True)
    (output / 'summary.json').write_text(json.dumps(report, indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), layout='constrained')
    names = ['1-15', '16-99', '100-999', '1000+']
    axes[0].bar(names, [100 * r['object_fraction'] for r in grouped])
    axes[0].set(title='Object snapshot count', xlabel='Area (grid cells)', ylabel='Percent of labeled objects')
    axes[1].bar(names, [100 * (r['labeled_rain_fraction'] or 0) for r in grouped])
    axes[1].set(title='Labeled rain contribution', xlabel='Area (grid cells)', ylabel='Percent of labeled mm-cells')
    axes[2].hist(spans, bins=np.arange(.5, min(48, max(spans, default=1)) + 1.5))
    axes[2].set(title='Observed branch spans (<=48 h shown)', xlabel='Hours, gaps included', ylabel='Branches', yscale='log')
    fig.savefig(output / 'rain_size_and_lifetimes.png', dpi=150)
    plt.close(fig)
    (output / 'REPORT.md').write_text('# Rain-catalog screening\n\n' + report['limitations'] +
        '\n\n![Size and lifetime screening](rain_size_and_lifetimes.png)\n\n```json\n' +
        json.dumps(report, indent=2) + '\n```\n')
    (output / 'SUCCESS').write_text('Catalog screening completed; not accuracy validation\n')
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path)
    p.add_argument('--output-dir', type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(summarize(a.source, a.output_dir), indent=2))
