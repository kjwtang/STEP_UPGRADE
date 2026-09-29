#!/usr/bin/env python3
"""Inspect saved streaming gap edges without rerunning identification/tracking.

Unadvected endpoint-footprint overlap is diagnostic evidence, not a replay of
the tracker's candidate gates or proof of meteorological continuity.
"""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def footprint_candidates(labels, footprint, lookup):
    values, counts = np.unique(labels[footprint], return_counts=True)
    result = []
    for label, count in zip(values, counts):
        if label == 0:
            continue
        node = lookup[int(label)]
        result.append(dict(node_id=int(node['node_id']),
                           local_label=int(label), area_cells=int(node['area_cells']),
                           intersection_cells=int(count)))
    return sorted(result, key=lambda item: (-item['intersection_cells'], item['node_id']))


def audit(source, output, plot=True):
    nodes, by_time, edges, frames = {}, defaultdict(dict), [], {}
    # Only chunk-local catalogs: never also read a finalized aggregate catalog.
    for part in sorted(source.glob('chunk_*')):
        if not part.is_dir():
            continue
        labels = np.load(part / 'identified_labels.npy', mmap_mode='r')
        offset = int(part.name.split('_')[-1])
        for i in range(len(labels)):
            if offset+i in frames:
                raise ValueError('Overlapping chunks')
            frames[offset+i] = labels[i]
        for node in rows(part / 'objects.csv'):
            nid, time = int(node['node_id']), int(node['time'])
            if nid in nodes:
                raise ValueError(f'Duplicate node {nid}')
            nodes[nid] = node
            by_time[time][int(node['local_label'])] = node
        edges.extend(rows(part / 'edges.csv'))
    if not frames:
        raise ValueError('No saved chunk labels found; run requires --save-labels')
    incoming, outgoing = defaultdict(list), defaultdict(list)
    for edge in edges:
        incoming[int(edge['child_node_id'])].append(edge)
        outgoing[int(edge['parent_node_id'])].append(edge)
    output.mkdir(parents=True, exist_ok=False)
    gaps = [e for e in edges if int(e['gap']) > 0]
    report = []
    for index, edge in enumerate(gaps, 1):
        parent, child = nodes[int(edge['parent_node_id'])], nodes[int(edge['child_node_id'])]
        t0, t1 = int(parent['time']), int(child['time'])
        if t1-t0 != int(edge['gap'])+1:
            raise ValueError('Inconsistent gap duration')
        pm = frames[t0] == int(parent['local_label'])
        cm = frames[t1] == int(child['local_label'])
        if pm.sum() != int(parent['area_cells']) or cm.sum() != int(child['area_cells']):
            raise ValueError('Endpoint labels and catalog areas disagree')
        record = dict(edge=edge, parent=parent, child=child, intermediate=[])
        for time in range(t0+1, t1):
            candidates = footprint_candidates(frames[time], pm | cm, by_time[time])
            for candidate in candidates:
                nid = candidate['node_id']
                candidate['incoming_edges'] = incoming[nid]
                candidate['outgoing_edges'] = outgoing[nid]
            record['intermediate'].append(dict(time=time, candidates=candidates))
        report.append(record)
        if plot:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            yy, xx = np.where(pm | cm)
            y0, y1 = max(0, yy.min()-15), min(pm.shape[0], yy.max()+16)
            x0, x1 = max(0, xx.min()-15), min(pm.shape[1], xx.max()+16)
            fig, axes = plt.subplots(1, t1-t0+1, figsize=(5*(t1-t0+1), 5), squeeze=False)
            for time, ax in zip(range(t0, t1+1), axes[0]):
                crop = frames[time][y0:y1, x0:x1]
                ax.imshow(crop > 0, origin='lower', cmap='Greys', vmin=0, vmax=1)
                # Fixed endpoint footprints over all panels; not advected.
                for mask, color in ((pm, 'red'), (cm, 'blue')):
                    local = mask[y0:y1, x0:x1]
                    if local.any() and not local.all():
                        ax.contour(local, levels=[.5], colors=[color], linewidths=.6)
                for label in np.unique(crop):
                    if label:
                        y, x = np.where(crop == label)
                        ax.text(x.mean(), y.mean(), str(label), fontsize=6)
                ax.set_title(f'Frame {time}; numbers = local labels')
            fig.suptitle(f"Gap {parent['node_id']} -> {child['node_id']} | red: parent; blue: child; black: identified rain")
            fig.tight_layout()
            fig.savefig(output / f'gap_{t0:03d}_{parent["node_id"]}_{child["node_id"]}.png', dpi=130)
            plt.close(fig)
        print(f'[{index}/{len(gaps)}] gap {t0}->{t1}', flush=True)
    (output / 'gap_audit.json').write_text(json.dumps(report, indent=2))
    (output / 'README.txt').write_text(
        'Diagnostic only: black is identified rain, not rain intensity. Red/blue are fixed endpoint footprints.\n'
        'Intermediate candidates overlap either unadvected endpoint footprint; absence does not prove no moving counterpart.\n'
        'Candidate edges describe saved graph decisions, not the rejected scores/gates. No parameters or algorithms changed.\n')
    print(f'Completed: {len(report)} gaps; {output}')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--no-plots', action='store_true')
    args = parser.parse_args()
    audit(args.source, args.output_dir, not args.no_plots)
