#!/usr/bin/env python3
"""Audit event screening and threshold sensitivity on saved identification.

Does not change production defaults or rerun identification/original tracking.
Candidate rows are streamed to disk; each threshold has independent causal state.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from step.tracking import Tracker, TrackGraph, _objects
from run_npy_validation import write_graph
from validate_real_data import provenance


def structure_counts(candidates):
    counts = Counter()
    for parents, children in Tracker._components(candidates):
        if len(parents) == len(children) == 1:
            continue
        counts['split' if len(parents) == 1 else
               'merge' if len(children) == 1 else 'complex'] += 1
    return dict(counts)


def audit(labels, rain, directory, tau, threshold, displacement, saved=None,
          event_policy='score_and_overlap', event_overlap=.10, event_min_pixels=1):
    directory.mkdir()
    tracker = Tracker(tau=tau, max_displacement=displacement, max_gap=1,
                      gap_tau=.45, gap_ambiguity=.05, event_overlap=event_overlap,
                      sequence_id='original_comparison', event_policy=event_policy,
                      event_min_pixels=event_min_pixels)
    combined = TrackGraph()
    frames = []
    fields = ['frame', 'parent_frame', 'parent_label', 'child_label', 'score',
              'parent_coverage', 'child_coverage', 'passes_tau',
              'passes_event_overlap', 'passes_event_policy', 'accepted_event']
    started = time.perf_counter()
    with (directory/'candidates.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for t in range(len(rain)):
            data = np.where(rain[t] >= threshold, rain[t], 0).astype(np.float32)
            parents = list(tracker.state.active)
            children = _objects(labels[t], data)
            candidates = tracker._candidates(parents, children, t)
            overlap = [c for c in candidates if
                       max(c.parent_coverage, c.child_coverage) >= tracker.event_overlap]
            eligible = tracker._event_candidates(candidates, parents, children)
            eligible_pairs = {(c.parent_index,c.child_index) for c in eligible}
            raster, graph = tracker.update(t, labels[t], data)
            if saved is not None and not np.array_equal(raster, saved[t]):
                raise ValueError(f'Baseline replay differs at frame {t}; use revision3 results')
            child_labels = {n.node_id: n.object.label for n in graph.objects}
            accepted = {(e.parent_id, child_labels[e.child_id]): e.event
                        for e in graph.edges}
            for c in candidates:
                parent, child = parents[c.parent_index], children[c.child_index]
                writer.writerow(dict(frame=t, parent_frame=parent.node.time,
                    parent_label=parent.node.object.label, child_label=child.label,
                    score=c.score, parent_coverage=c.parent_coverage,
                    child_coverage=c.child_coverage, passes_tau=c.score >= tau,
                    passes_event_overlap=max(c.parent_coverage,c.child_coverage) >= event_overlap,
                    passes_event_policy=(c.parent_index,c.child_index) in eligible_pairs,
                    accepted_event=accepted.get((parent.node.node_id,child.label),'')))
            frames.append(dict(frame=t, candidates=len(candidates),
                overlap_structures_before_tau=structure_counts(overlap),
                event_structures_after_tau=structure_counts([c for c in overlap if c.score >= tau]),
                event_structures_after_policy=structure_counts(eligible),
                accepted_edge_types=dict(Counter(e.event for e in graph.edges))))
            combined.objects.extend(graph.objects)
            combined.edges.extend(graph.edges)
            combined.events.extend(graph.events)
            combined.family_map = graph.family_map
    write_graph(combined, directory)
    result = dict(tau=tau, gap_tau=.45, event_overlap=event_overlap,
        directory=str(directory.resolve()),
        event_policy=event_policy, event_min_pixels=event_min_pixels,
        baseline_raster_equal=True if saved is not None else None,
        audit_wall_seconds=time.perf_counter()-started,
        edge_types=dict(Counter(e.event for e in combined.edges)),
        event_types=dict(Counter(e.event for e in combined.events)), frames=frames)
    (directory/'summary.json').write_text(json.dumps(result, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('comparison', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--taus', nargs='+', type=float, default=[.30,.325,.35])
    parser.add_argument('--include-overlap-policy', action='store_true',
                        help='Also test overlap-only events; baseline is always verified first')
    parser.add_argument('--overlap-policy-fraction', type=float, default=.5)
    parser.add_argument('--overlap-policy-min-pixels', type=int, default=2)
    args = parser.parse_args()
    if any(not 0 <= x <= 1 for x in args.taus):
        parser.error('taus must be between 0 and 1')
    if not 0 < args.overlap_policy_fraction <= 1 or args.overlap_policy_min_pixels < 1:
        parser.error('Invalid overlap policy settings')
    source = args.comparison.resolve()
    config = json.loads((source/'configuration.json').read_text())
    spec = json.loads((source/'new_track/spec.json').read_text())
    rain = np.load(source/'rain_mm_h.npy', mmap_mode='r')
    labels = np.load(source/'new_id/labels.npy', mmap_mode='r')
    saved = np.load(source/'new_track/labels.npy', mmap_mode='r')
    if rain.shape != labels.shape or labels.shape != saved.shape:
        raise ValueError('Saved shapes differ')
    digest = hashlib.sha256()
    for frame in rain:
        digest.update(frame.tobytes())
    if digest.hexdigest() != config['input_sha256']:
        raise ValueError('Saved rain hash differs')
    args.output_dir.mkdir(parents=True, exist_ok=False)
    (args.output_dir/'provenance.json').write_text(json.dumps(provenance(args),indent=2,default=str))
    # Verify baseline first, before interpreting sensitivity runs.
    baseline = spec['new_tau']
    taus = [baseline] + [x for x in sorted(set(args.taus)) if x != baseline]
    results = []
    for tau in taus:
        print(f'Auditing all frames: tau={tau}', flush=True)
        results.append(audit(labels,rain,args.output_dir/f'tau_{tau}',tau,
            config['threshold'],spec['new_displacement'],saved if tau==baseline else None))
    if args.include_overlap_policy:
        print('Auditing matched-overlap score control',flush=True)
        results.append(audit(labels,rain,args.output_dir/'score_policy_matched_overlap',baseline,
            config['threshold'],spec['new_displacement'],
            event_overlap=args.overlap_policy_fraction,
            event_min_pixels=args.overlap_policy_min_pixels))
        print('Auditing experimental overlap-only events',flush=True)
        results.append(audit(labels,rain,args.output_dir/'overlap_policy',baseline,
            config['threshold'],spec['new_displacement'],event_policy='overlap',
            event_overlap=args.overlap_policy_fraction,
            event_min_pixels=args.overlap_policy_min_pixels))
    (args.output_dir/'summary.json').write_text(json.dumps(results,indent=2))
    from recheck_tracking import edge_map
    baseline_edges=edge_map(Path(results[0]['directory']))
    changes=[]
    for result in results[1:]:
        edges=edge_map(Path(result['directory']))
        changes.append(dict(directory=result['directory'],
            added=[dict(object_pair=k,event=edges[k]) for k in sorted(edges.keys()-baseline_edges.keys())],
            removed=[dict(object_pair=k,event=baseline_edges[k]) for k in sorted(baseline_edges.keys()-edges.keys())],
            changed_event=[dict(object_pair=k,old=baseline_edges[k],new=edges[k])
                for k in sorted(edges.keys()&baseline_edges.keys()) if edges[k]!=baseline_edges[k]]))
    (args.output_dir/'edge_changes_vs_baseline.json').write_text(json.dumps(changes,indent=2))
    (args.output_dir/'README.txt').write_text(
        'Overlap structures before tau are hypotheses, not true events. Candidate CSVs '
        'contain adjacent active pairs only; gap edges are in edges.csv. Gap threshold '
        'remains .45 throughout, so this is not a gap-threshold sweep. Compare links '
        'by frame/local-label through objects.csv, not branch ID. Audit timing includes '
        'duplicate candidate evaluation and CSV I/O; it is not a tracking benchmark. '
        'All-continue output alone does not prove an error. Review spatial evidence '
        'before choosing thresholds. Production defaults are unchanged.\n')
    (args.output_dir/'SUCCESS').touch()
    print(json.dumps([{k:v for k,v in r.items() if k!='frames'} for r in results],indent=2))


if __name__ == '__main__':
    main()
