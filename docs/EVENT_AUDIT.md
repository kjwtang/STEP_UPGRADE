# Event audit and hierarchical connectivity

Revision 3 RCC replay recovered four continuation edges (140 to 144), with no
removed edges. All edges were continuations. This does not establish that real
splits, merges or gaps exist in those 24 hours, nor that event detection is broken.
Events are screened before one-to-one assignment, not after it.

Run on the completed revision 3 output, in the existing compute-node environment:

```bash
python -u scripts/audit_tracking_events.py \
  results/tracking_revision3_600 \
  --output-dir results/event_audit_600
```

This reuses precipitation and identified objects, verifies the rain hash and
baseline raster replay, and runs independent causal trackers at tau .35, .30,
and .325. The baseline is checked first. No production defaults change.
Each run exports objects, edges, events, family map, a streamed candidates CSV,
and frame-level summaries. Review overlap-supported nontrivial components before
and after score filtering: these indicate *hypothetical* split/merge structures,
not meteorological truth. All candidate rows describe adjacent active objects;
gap edges are exported but dormant-candidate rejection reasons are not audited.
Gap tau remains .45; this experiment does not sweep gap parameters.

Compare accepted edges by parent/child frame and local identification label using
objects.csv, not branch numbers. Inspect the spatial masks for newly accepted
links and for components lost at the score threshold. Lower tau can also change
event topology and subsequent velocities, so outcomes need not be nested or
monotonic. Audit runtime includes extra scoring and CSV writes and is not the
production tracking runtime. SUCCESS is written only after all runs finish.

## What to borrow from HDBSCAN

HDBSCAN constructs a hierarchy using mutual-reachability distances and a minimum
spanning tree, then selects stable clusters:
https://hdbscan.readthedocs.io/en/latest/how_hdbscan_works.html

The useful analogy is examining connectivity across thresholds rather than
tuning one threshold to one storm. Our three-run sweep is only a sensitivity
experiment, NOT HDBSCAN, its formal cluster stability, or a calibrated confidence.

For identification, a separate experimental hierarchy over rain objects/cores
could represent a broad system containing several cores; physical distances,
rain support and weak bridges would need validation against the frozen baseline.
For tracking, keep the time-directed graph: a static undirected clustering tree
does not itself encode continuation, split then merge, or missing frames. Do not
replace the graph by an MST, which preserves threshold connectivity but not all
lineage edges. Do not flatten all high-resolution hourly wet pixels into one
unbounded clustering problem. No HDBSCAN dependency or algorithm is added here.
