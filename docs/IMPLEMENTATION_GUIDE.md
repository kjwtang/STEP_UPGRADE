# STEP_UPGRADE 0.2 implementation and production guide

This document records the behavior that downstream implementations and RCC
production wrappers must preserve. It separates the implemented 0.2 core from
work that still requires real-data calibration.

## Scope

The implemented core provides:

1. deterministic binary-mask identification with thresholds and missing-data
   barriers;
2. unique object nodes and uninterrupted branch IDs;
3. continuation, split, merge and complex event edges;
4. one-frame or configurable below-threshold gap reconnection;
5. portable JSON state for continuous chunk/month processing;
6. node, edge, event and branch-to-family outputs.

CAPE, updraft, ECAPE and microphysical sampling belong in a later module that
uses node masks and branch/family IDs. Do not add those fields to the matching
score without a separate scientific decision and validation.

## Input contract

- Precipitation is `(time, y, x)` and must have a documented interval and unit.
- Existing monthly `PRECIP(t)` data represent `(t-1h,t]` in mm and must not be
  differenced again.
- Raw WRF `RAINNC` is cumulative and needs a separate, restart-aware conversion
  before use.
- NaN or an explicit false `valid_mask` value is missing, not zero rain.
- Grid, projection, coordinates and cadence must be constant within a sequence.
- Present 2005 and future 2059 are independent sequences and never share state.
- A month or file boundary is not a storm boundary.

Before production, build an ordered manifest and reject duplicate, decreasing,
or unexpectedly missing timestamps. If a whole observation frame is missing,
the current public API needs the wrapper to decide whether to stop or advance
with an explicit empty but valid frame. Do not silently substitute a missing
file with zero rain.

## Identification behavior

`identify()` computes a rain mask using finite/valid data, an explicit
threshold, and positive precipitation. It applies binary dilation, propagates
only through valid pixels, performs 8-connected labeling, and projects labels
back onto original rain pixels. `min_size` counts original rain pixels, not the
dilated bridge.

Labels are sorted by first row-major rain pixel. Parallel worker count and
chunk size therefore cannot change them. The morphology parameter is a
dilation radius/footprint around each rain pixel. Two objects can join over a
separation approaching twice the radius; do not describe radius as the maximum
inter-object gap.

The historical production baseline was threshold 1 mm h-1, radius 9 cells and
min_size 1. Treat it as a reproducibility point rather than a calibrated final
choice. Evaluate at least threshold `[0.6, 1, 2]` mm h-1 and radius
`[0, 1, 3, 5, 9]` cells on representative cases, holding other parameters
fixed while interpreting each change.

## Lineage semantics

Each observed two-dimensional object gets one immutable `node_id`. A
`branch_id` continues only through accepted one-to-one continuation or
one-to-one gap continuation. At split, merge or complex events all incoming
branches end and every child receives a new branch. This prevents disconnected
children from sharing one persistent raster ID.

All branches connected by accepted event edges form a family. The state stores
a union-find mapping. A future merge can combine old families, so consumers
must use the final canonical `family_map` rather than an early node's
provisional family value.

Matching first predicts the next centroid from the previous two observations.
A KD-tree limits candidates to the configured displacement per elapsed frame.
For every candidate the tracker saves raw and motion-advected IoU, predicted
distance, intensity continuity and directional coverage. Event components need
overlap evidence. Non-event candidates use deterministic maximum-score
one-to-one assignment.

The current `km` parameter is retained for compatibility but means grid cells
per frame. Production wrappers must convert a physical speed and timestep into
this value and record the conversion. `phi` is accepted but ignored.

## Gap behavior

`max_gap=1` allows an object last observed at `t` to connect to an object at
`t+2`; it does not permit two missing object frames. Gap matching is
one-to-one, uses `gap_tau`, and never creates an interpolated node or rainfall.
Branches that participate in split/merge/complex events are consumed and are
not placed into the dormant pool.

The current implementation models below-threshold object gaps. Missing or
invalid full frames need explicit wrapper policy and should be recorded
separately before scientific production. Domain-edge entry/exit censoring and
gap ambiguity rejection are recommended follow-up metadata if the real-data
review shows they materially affect results.

## Continuous chunk and month processing

Tracking is chronological and serial within one sequence; identification can
be parallel. Call `track_with_graph(..., return_state=True)`, atomically save
the returned state, and pass it into the next chunk or month. Do not invoke
tracking separately for each Slurm array month.

State includes active/dormant node masks, branch and event counters, motion
history, last time, and family union state. `save_tracking_state()` writes a
temporary JSON file, flushes it, and atomically replaces the target. A
production wrapper should similarly write chunk outputs into a temporary
generation and publish its manifest only after tables, rasters, checksums and
the matching checkpoint are complete.

Completed graph nodes retain statistics but discard pixel coordinate lists;
only active/dormant state keeps the masks needed for the next match. This keeps
memory proportional to the current tracking frontier rather than the full
season's rain-pixel history.

Recommended RCC scheduling choices are one sequential multi-month job or a
month chain with `afterok`, where each job reads the preceding state. Different
climates, ensemble members, or parameter experiments may run independently.

## Output tables

The validation runner produces:

- `identified_labels.npy`: local object labels per frame;
- `tracked_labels.npy`: branch IDs per rain pixel;
- `objects.csv`: node, branch, family and object statistics;
- `edges.csv`: accepted lineage evidence and scores;
- `events.csv`: split, merge and complex event memberships;
- `family_map.csv`: current canonical branch-to-family mapping;
- `run_settings.json`: exact invocation and shapes;
- optional JSON state: continuation input for the next chunk/month.

For a final season catalog, concatenate node/edge/event partitions, then apply
the final family map to every branch. Do not sum rainfall along all family
paths without defining how split/merge conservation is handled; node-level
precipitation must be counted once.

`scripts/finalize_catalog.py` performs this table merge, validates unique node
IDs and edge references, and writes the canonical family IDs. Supply every
chronological part of the sequence; omitting an earlier part makes cross-part
edge validation fail.

## Required tests and review

The checked-in tests must pass. Extend them when behavior changes. Required
invariants are:

- edge times increase and no graph cycles occur;
- one branch appears at most once per frame;
- split/merge children use new branches;
- each node/local-label mapping is unique;
- accepted edges reference existing nodes;
- a terminal and child participate in at most one gap continuation;
- whole-run and chunk/checkpoint runs yield identical raster IDs, nodes, edges
  and events.

Before full JJA processing, inspect multiple 24–72 hour animations containing
isolated convection, organized systems, weak/stratiform rain, fast motion,
month boundaries, and both climates. Compare old and new identification,
review every split/merge in selected cases, and report object counts, lifetime
distribution, branch speeds, event frequency, gap frequency, largest-family
size, runtime and memory. The legacy tracker is a comparison, not truth.

Production may start only after the real-data review selects threshold,
dilation, displacement, score, overlap and gap settings. Preserve the tested
configuration and repository commit in every output manifest.
