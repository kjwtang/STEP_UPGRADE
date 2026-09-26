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


def audit(labels, rain, directory, tau, threshold, displacement, saved=None):
    directory.mkdir()
    tracker = Tracker(tau=tau, max_displacement=displacement, max_gap=1,
                      gap_tau=.45, gap_ambiguity=.05, event_overlap=.10,
                      sequence_id='original_comparison')
    combined = TrackGraph()
    frames = []
    fields = ['frame', 'parent_frame', 'parent_label', 'child_label', 'score',
              'parent_coverage', 'child_coverage', 'passes_tau',
              'passes_event_overlap', 'accepted_event']
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
            eligible = [c for c in overlap if c.score >= tau]
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
                    passes_event_overlap=max(c.parent_coverage,c.child_coverage) >= .10,
                    accepted_event=accepted.get((parent.node.node_id,child.label),'')))
            frames.append(dict(frame=t, candidates=len(candidates),
                overlap_structures_before_tau=structure_counts(overlap),
                event_structures_after_tau=structure_counts(eligible),
                accepted_edge_types=dict(Counter(e.event for e in graph.edges))))
            combined.objects.extend(graph.objects)
            combined.edges.extend(graph.edges)
            combined.events.extend(graph.events)
            combined.family_map = graph.family_map
    write_graph(combined, directory)
    result = dict(tau=tau, gap_tau=.45, event_overlap=.10,
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
    args = parser.parse_args()
    if any(not 0 <= x <= 1 for x in args.taus):
        parser.error('taus must be between 0 and 1')
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
    (args.output_dir/'summary.json').write_text(json.dumps(results,indent=2))
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
