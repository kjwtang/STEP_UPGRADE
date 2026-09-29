#!/usr/bin/env python3
"""Replay saved object statistics/masks and explain selected adjacent links.

No NetCDF or identification rerun. This validates tracking replay, NOT original
rain ingestion or the saved object statistics themselves.
"""
import argparse
import csv
from dataclasses import asdict
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
from step.tracking import Tracker, StormObject, TrackGraph, load_tracking_state
from diagnose_tracking import evidence
from run_npy_validation import write_graph

DEFAULT_PAIRS = [(7414,7500),(7500,7582),(7582,7669),(7669,7761),
                 (7624,7722),(7722,7807),(7807,7927)]


def read_rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def restore_objects(labels, rows):
    result = []
    if set(np.unique(labels)) - {0} != {int(r['local_label']) for r in rows}:
        raise ValueError('Label/catalog disagreement')
    for r in sorted(rows, key=lambda r: int(r['local_label'])):
        label = int(r['local_label'])
        y0,y1,x0,x1 = (int(r[k]) for k in ('y_start','y_stop','x_start','x_stop'))
        yy,xx = np.where(labels[y0:y1,x0:x1] == label)
        if len(yy) != int(r['area_cells']):
            raise ValueError('Area mismatch')
        result.append(StormObject(label, len(yy),
            (float(r['centroid_y']),float(r['centroid_x'])), float(r['mean_intensity']),
            ((y0,y1),(x0,x1)), float(r['max_intensity']), float(r['precipitation_sum']),
            list(zip((yy+y0).tolist(),(xx+x0).tolist()))))
    return result


def replay(source, output, pairs):
    if output.exists():
        raise FileExistsError(output)
    parts = sorted(source.glob('chunk_*'))
    if not parts:
        raise ValueError('No chunks')
    final = load_tracking_state(parts[-1]/'state.json')
    c = final.tracking_config
    tracker = Tracker(tau=c['tau'], max_displacement=c['max_displacement'],
        max_gap=int(c['max_gap']), gap_tau=c['gap_tau'], gap_ambiguity=c['gap_ambiguity'],
        event_overlap=c['event_overlap'], sequence_id=final.sequence_id,
        event_policy='overlap' if c.get('event_policy_overlap') else 'score_and_overlap',
        event_min_pixels=int(c.get('event_min_pixels',1)),
        gap_conflict_policy='endpoint_overlap' if c.get('gap_conflict_endpoint_overlap_v1') else 'off')
    if tracker.state.tracking_config != c:
        raise ValueError('Unsupported saved algorithm/configuration')
    catalog = {int(r['node_id']):r for part in parts for r in read_rows(part/'objects.csv')}
    for p,ch in pairs:
        if p not in catalog or ch not in catalog:
            raise ValueError(f'Missing requested pair {p}->{ch}')
        if int(catalog[ch]['time']) != int(catalog[p]['time'])+1:
            raise ValueError('Pairs must be adjacent')
    requested = {}
    for pair in pairs:
        requested.setdefault(int(catalog[pair[1]]['time']), []).append(pair)
    decisions, conflicts, checked = [], [], 0
    for part in parts:
        labels = np.load(part/'identified_labels.npy', mmap_mode='r')
        saved = np.load(part/'tracked_labels.npy', mmap_mode='r')
        offset = int(part.name.split('_')[-1])
        combined = TrackGraph()
        for i, frame in enumerate(labels):
            t = offset+i
            if t != checked:
                raise ValueError('Chunks must start at zero and be contiguous')
            rows = [r for r in catalog.values() if int(r['time']) == t]
            children = restore_objects(frame, rows)
            parents = list(tracker.state.active)
            pending = []
            if t in requested:
                candidates = tracker._candidates(parents, children, t)
                cmap = {(x.parent_index,x.child_index):x for x in candidates}
                event_pairs = {(x.parent_index,x.child_index) for x in
                               tracker._event_candidates(candidates, parents, children)}
                for p,ch in requested[t]:
                    pi = next(j for j,v in enumerate(parents) if v.node.node_id == p)
                    ci = next(j for j,v in enumerate(children) if v.label == int(catalog[ch]['local_label']))
                    e = evidence(tracker, parents[pi], children[ci], t)
                    candidate = cmap.get((pi,ci))
                    if candidate and not np.isclose(candidate.score,e['score_if_evaluated'],rtol=0,atol=1e-12):
                        raise ValueError('Score formula no longer reproduces implementation')
                    pending.append(dict(parent_node_id=p,child_node_id=ch,time=t,
                        admitted=candidate is not None,event_candidate=(pi,ci) in event_pairs,
                        area_ratio=children[ci].area/parents[pi].node.object.area,**e))
            # _objects is replaced only for this isolated replay call. All
            # tracker decisions still use its actual production implementation.
            with patch('step.tracking._objects', return_value=children):
                raster, graph = tracker.update(t, frame, np.zeros(frame.shape, dtype=np.float32))
            if not np.array_equal(raster,saved[i]):
                raise ValueError(f'Raster replay mismatch at {t}; do not interpret diagnostics')
            for row in pending:
                p,ch = row['parent_node_id'],row['child_node_id']
                accepted = [e for e in graph.edges if e.parent_id==p and e.child_id==ch]
                competitors = [asdict(e) for e in graph.edges if e.parent_id==p or e.child_id==ch]
                row['accepted_edges'] = [asdict(e) for e in accepted]
                row['endpoint_edges'] = competitors
                row['decision'] = ('accepted' if accepted else 'outside_candidate_gates' if not row['admitted']
                    else 'below_score_threshold' if not row['passes_score']
                    else 'endpoint_used_by_other_edge' if competitors else 'eligible_unassigned')
                decisions.append(row)
            conflicts.extend(dict(**r, passes_gap_score=r['score']>=tracker.gap_tau)
                             for r in tracker.gap_conflicts)
            combined.objects.extend(graph.objects)
            combined.edges.extend(graph.edges)
            combined.events.extend(graph.events)
            combined.family_map = graph.family_map
            checked += 1
            print(f'[{checked}/{final.last_time+1}] replay verified', flush=True)
        with TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            write_graph(combined,tmp)
            for name in ('objects.csv','edges.csv','events.csv','family_map.csv'):
                if (tmp/name).read_bytes() != (part/name).read_bytes():
                    raise ValueError(f'Graph replay mismatch: {part.name}/{name}')
    report = dict(replay_equal=True,frames_replayed=checked,tracking_config=c,
                  pairs=decisions,gap_conflicts=conflicts,
                  note='Saved statistics/mask replay, not independent ingestion validation. Outside-gate scores are counterfactual; reasons are pair-level, not causal attribution.')
    output.mkdir(parents=True)
    (output/'diagnosis.json').write_text(json.dumps(report,indent=2))
    lines = ['# Adjacent-link replay diagnosis','',f'Exact raster and catalog replay: {checked} frames.','']
    for r in decisions:
        lines.append(f"- {r['parent_node_id']} -> {r['child_node_id']}: {r['decision']}; "
                     f"score={r['score_if_evaluated']:.4f}; raw/adv IoU={r['raw_iou']:.4f}/{r['advected_iou']:.4f}; "
                     f"distance/prediction error={r['actual_displacement_cells']:.2f}/{r['prediction_error_cells']:.2f}; "
                     f"area ratio={r['area_ratio']:.2f}; event candidate={r['event_candidate']}")
    lines += ['',report['note']]
    (output/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path)
    parser.add_argument('--output-dir',required=True,type=Path)
    parser.add_argument('--pair',action='append',help='Parent:child node IDs; repeatable')
    args = parser.parse_args()
    pairs = [tuple(map(int,p.split(':'))) for p in args.pair] if args.pair else DEFAULT_PAIRS
    replay(args.source,args.output_dir,pairs)
