#!/usr/bin/env python3
"""Replay saved object statistics/masks and explain selected adjacent links.

No NetCDF or identification rerun. This validates tracking replay, NOT original
rain ingestion or the saved object statistics themselves.
"""
import argparse
import csv
from contextlib import contextmanager
from dataclasses import asdict
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import threading
import time
from unittest.mock import patch

import numpy as np
from step.tracking import Tracker, StormObject, TrackGraph, load_tracking_state, save_tracking_state, _pixel_overlap
from diagnose_tracking import evidence
from run_npy_validation import write_graph

DEFAULT_PAIRS = [(7414,7500),(7500,7582),(7582,7669),(7669,7761),
                 (7624,7722),(7722,7807),(7807,7927)]


@contextmanager
def timed_stage(label, timings):
    """Heartbeat exposes long frame-internal work; it is not an ETA."""
    start=time.perf_counter()
    stop=threading.Event()
    print(f'{label}: START',flush=True)
    def heartbeat():
        while not stop.wait(15):
            print(f'{label}: still running ({time.perf_counter()-start:.1f}s)',flush=True)
    thread=threading.Thread(target=heartbeat,daemon=True)
    thread.start()
    ok=False
    try:
        yield
        ok=True
    finally:
        stop.set()
        thread.join()
        elapsed=time.perf_counter()-start
        timings.append(dict(stage=label,seconds=elapsed,completed=ok))
        print(f'{label}: {"DONE" if ok else "INTERRUPTED/FAILED"} {elapsed:.3f}s',flush=True)


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


def replay(source, output, pairs, checkpoint_every=0):
    if checkpoint_every < 0:
        raise ValueError('checkpoint_every must be nonnegative')
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
        gap_conflict_policy='endpoint_overlap' if c.get('gap_conflict_endpoint_overlap_v1') else 'off',
        adjacent_policy='overlap_first' if c.get('adjacent_overlap_first_v1') else 'score',
        velocity_reset=('morphology' if c['velocity_reset_morphology_v1'] else 'off')
                       if 'velocity_reset_morphology_v1' in c else 'auto',
        score_policy='coherent_adjacent' if c.get('coherent_adjacent_score_v1') else 'legacy')
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
    timings=[]
    by_time={}
    for r in catalog.values():
        by_time.setdefault(int(r['time']),[]).append(r)
    for part in parts:
        labels = np.load(part/'identified_labels.npy', mmap_mode='r')
        saved = np.load(part/'tracked_labels.npy', mmap_mode='r')
        offset = int(part.name.split('_')[-1])
        combined = TrackGraph()
        for i, frame in enumerate(labels):
            t = offset+i
            if t != checked:
                raise ValueError('Chunks must start at zero and be contiguous')
            label=f'frame {t+1}/{final.last_time+1} ({100*(t+1)/(final.last_time+1):.1f}%)'
            with timed_stage(label+' restore objects',timings):
                children = restore_objects(frame, by_time.get(t,[]))
            parents = list(tracker.state.active)
            pending = []
            if t in requested:
                for p,ch in requested[t]:
                    pi = next(j for j,v in enumerate(parents) if v.node.node_id == p)
                    ci = next(j for j,v in enumerate(children) if v.label == int(catalog[ch]['local_label']))
                    # Admission and pair eligibility are pair-local. Evaluate
                    # only this pair; full competition still runs in update.
                    with timed_stage(label+f' pair {p}->{ch}',timings):
                        candidates=tracker._candidates([parents[pi]],[children[ci]],t)
                        strong=tracker._strong_overlap_candidates(candidates,[parents[pi]],[children[ci]])
                        event=tracker._event_candidates(candidates,[parents[pi]],[children[ci]])
                        terminal = parents[pi]
                        raw = _pixel_overlap(terminal.node.object, children[ci])
                        velocity = tracker._velocity(terminal)
                        shift = tuple(int(round(v*(t-terminal.node.time))) for v in velocity)
                        adv = _pixel_overlap(terminal.node.object, children[ci], shift)
                        _,rpc,rcc=raw
                        _,apc,acc=adv
                        e = evidence(tracker, terminal, children[ci], t,
                                     raw_overlap=raw,advected_overlap=adv)
                    candidate = candidates[0] if candidates else None
                    if candidate and not np.isclose(candidate.score,e['score_if_evaluated'],rtol=0,atol=1e-12):
                        raise ValueError('Score formula no longer reproduces implementation')
                    pending.append(dict(parent_node_id=p,child_node_id=ch,time=t,
                        admitted=candidate is not None,event_candidate=bool(event or strong),
                        strong_raw_overlap=bool(strong),
                        object_pair=[int(catalog[p]['time']),int(catalog[p]['local_label']),
                                     int(catalog[ch]['time']),int(catalog[ch]['local_label'])],
                        raw_parent_coverage=rpc,raw_child_coverage=rcc,
                        advected_parent_coverage=apc,advected_child_coverage=acc,
                        raw_intersection_cells=round(rpc*terminal.node.object.area),
                        advected_intersection_cells=round(apc*terminal.node.object.area),
                        previous_centroid=terminal.previous_centroid,previous_time=terminal.previous_time,
                        parent_centroid=terminal.node.object.centroid,child_centroid=children[ci].centroid,
                        area_ratio=children[ci].area/parents[pi].node.object.area,**e))
            # _objects is replaced only for this isolated replay call. All
            # tracker decisions still use its actual production implementation.
            with timed_stage(label+' tracking',timings):
                with patch('step.tracking._objects', return_value=children):
                    raster, graph = tracker.update(t, frame, np.zeros(frame.shape, dtype=np.float32))
            with timed_stage(label+' raster verification',timings):
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
            if checkpoint_every and checked % checkpoint_every == 0:
                with timed_stage(label+' checkpoint roundtrip',timings), TemporaryDirectory() as tmp:
                    path = Path(tmp)/'state.json'
                    save_tracking_state(tracker.state, path)
                    tracker.state = load_tracking_state(path)
            print(f'[{checked}/{final.last_time+1}] replay verified', flush=True)
        with timed_stage(part.name+' catalog verification',timings), TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            write_graph(combined,tmp)
            for name in ('objects.csv','edges.csv','events.csv','family_map.csv'):
                if (tmp/name).read_bytes() != (part/name).read_bytes():
                    raise ValueError(f'Graph replay mismatch: {part.name}/{name}')
    report = dict(replay_equal=True,frames_replayed=checked,tracking_config=c,
                  checkpoint_every=checkpoint_every,
                  pairs=decisions,gap_conflicts=conflicts,stage_timings=timings,
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
    parser.add_argument('--no-focus',action='store_true',help='Replay all frames without historical case-specific focus pairs')
    parser.add_argument('--checkpoint-every',type=int,default=0,help='Roundtrip state every N frames to verify another checkpoint boundary schedule')
    args = parser.parse_args()
    pairs = [tuple(map(int,p.split(':'))) for p in args.pair] if args.pair else ([] if args.no_focus else DEFAULT_PAIRS)
    replay(args.source,args.output_dir,pairs,args.checkpoint_every)
