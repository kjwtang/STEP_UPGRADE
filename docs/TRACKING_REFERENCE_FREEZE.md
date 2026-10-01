# Operational tracking reference freeze v1

Decision: freeze the **implementation and named reference configuration** for
subsequent identification/envelope experiments. Do not freeze a claim of
meteorological truth, original-STEP equivalence, convective classification, or
validated ten-year historical/future production. Scientific library defaults
remain unchanged; activate this reference explicitly.

Core commit: `daf70111e76fecfd90db69afb1aeeb8634526836`.
Preset: [rainfall_lineage_4km_1h_reference_v1.json](../configs/rainfall_lineage_4km_1h_reference_v1.json).
Its `tracking` dictionary can be passed to `Tracker` or `track_with_graph`;
identification and input-contract fields are separate, not tracker arguments.

## Frozen order

1. Significant raw split/merge/complex topology (balanced-v2 definition).
2. Strong symmetric raw-overlap continuations.
3. Score-qualified one-to-one assignment, with infeasible edges excluded first.
4. Mutual-dominant raw-overlap rescue **only for unmatched endpoints**.
5. Guarded gap recovery, ambiguity rejection, and 16-cell gap endpoint floor.

Ranked rescue uses >=16 intersection cells, >=15% raw coverage of both endpoints,
mutual best intersection among admitted candidates, and >=2x runner-up support.
Ranking includes reserved candidates; it does not fabricate dominance by removing
competitors. The rule is relational matching, not HDBSCAN. Constants are explicit
reference definitions, not calibrated physical probabilities. The displacement
parameter controls centroid search/scoring, not an absolute speed cap: overlapping
large masks can be admitted beyond their centroid radius.

Keep 1-mm/h identification and the 9-cell bridge fixed for this reference. Output
masks retain observed wet pixels only; morphology changes grouping. On a 4-km
grid a radius-9 dilation extends connectivity 36 km from each mask, so two nearby
pieces can bridge a substantially wider gap. This can chain rain bands; it is
not proof that one labeled object is one convective storm. Changing to a 0.1-mm/h
extent must be a separately versioned identification/rain-accounting experiment.

## Acceptance evidence

Full 1749 x 2049 grid, native hourly RAINNC differences. Six fixed windows below
total 528 processed intervals, with 48 overlapping June-30/July-1 intervals:
**480 distinct hours**. July-20 and August-25 were untouched holdouts when the
ranked-v2 rule was fixed. Other windows are development/regression, not holdouts.

| Window | Hours | Balanced-v2 edges | Frozen-reference edges | Added / removed | Reference gaps |
|---|---:|---:|---:|---:|---:|
| June 1 | 72 | 5,272 | 5,526 | 255 / 1 | 5 |
| June 30 | 72 | 4,714 | 5,074 | 361 / 1 | 16 |
| August 20 | 72 | 4,944 | 5,200 | 258 / 2 | 9 |
| July 1 | 168 | 11,809 | 12,522 | 717 / 4 | 38 |
| July 20 | 72 | 5,126 | 5,458 | 333 / 1 | 8 |
| August 25 | 72 | 4,696 | 4,948 | 255 / 3 | 12 |

All masks and object measurements match the balanced reference. No existing
continuation or event edge was removed or retyped in these comparisons. All 12
removed edges were gaps: three have replacement two-edge lineage paths, nine
do not. Other endpoint associations can change after rescue; fewer gaps are not
automatically better. Seven of eight preselected moderate-overlap interruptions
were recovered; one remained rejected. That small selected set is not an accuracy
estimate. Every large-object link with both endpoints >=100 cells has raw or
translated overlap in these windows; nearby small-object zero-overlap links
remain allowed by the scoring rule.

The July week reproduces all 168 tracked rasters and CSV catalogs with checkpoint
restoration every seven frames versus six. Structural audits pass unique nodes,
forward edges, timestamp/gap consistency, and connected, non-reused branch IDs.
107 tests pass, one is skipped. Synthetic tests include multi-way events,
ambiguity, expiry, ID stability and protection of score-qualified links from
weak rescue. Local NetCDF/NumPy environment warnings remain documented in
[the earlier assessment](LOCAL_FREEZE_ASSESSMENT.md).

Review masks and modified associations independently before production. These
tests establish reproducibility and targeted behavior, not observational truth.
Stationary-overlap rescue cannot solve every fast-moving small object's initial
association, nor every ambiguous group split. Core/updraft/CAPE classification
is outside the rain-only tracker.

## Worker evidence and ten-year execution design

Final reference, same 24 prepared full-grid frames, same process/environment:

| Frame workers | ID seconds | Sequential tracking seconds | ID + tracking |
|---:|---:|---:|---:|
| 1 | 10.84 | 3.59 | 14.43 |
| 4 | 3.12 | 3.55 | 6.67 |
| 8 | 1.75 | 3.45 | 5.20 |

Tracked raster SHA-256, full graph, scores, events, measurements and IDs match
across workers. Approximately 2.8x batch speedup at eight workers excludes read,
output and plotting, and is not an RCC forecast. The earlier two-round benchmark
also verifies halo tiles; frame parallelism was faster when enough frames were
available. Run the final reference benchmark independently on RCC:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -u scripts/benchmark_worker_modes.py PREPARED_MM_H.npy \
  --frames 32 --workers 1 4 8 16 --modes frames --repeat 2 \
  --adjacent-policy overlap_ranked --output-dir results/reference_workers_forward
python -u scripts/benchmark_worker_modes.py PREPARED_MM_H.npy \
  --frames 32 --workers 1 4 8 16 --modes frames --repeat 2 --reverse \
  --adjacent-policy overlap_ranked --output-dir results/reference_workers_reverse
```

For 32 workers, use at least 32 and preferably 64 frames per ID batch. Start with
16 CPUs / 32 GiB for bounded rain-only benchmarking; verify job-level memory
accounting, especially when increasing chunk length. Local child-memory sampling
lacked permissions, so parent RSS is not total node memory. 180 GiB is not a
requirement established by these tests.

Ten-year architecture: bounded rain ingestion -> parallel independent-frame ID
-> one ordered tracking-state stream -> chunk outputs and restart checkpoints.
Prioritize overlap of I/O and computation only with safe process isolation and
bounded queues; no asynchronous reader pipeline is claimed implemented here.
Do not fork a process while an I/O thread holds NetCDF/native-library state.
Do not independently track months/spatial tiles and concatenate IDs as if
equivalent. Parallelize truly independent climate/member sequences; hand state
forward where temporal continuity is required. A future exact boundary-stitcher
would require a proved algorithm, not merely overlap and ID renumbering.

Keep input/output storage out of Git. For production, avoid `--save-labels` on
every diagnostic run and plan checkpoint retention separately: current validation
saves every chunk for auditing, not a space-optimized multi-year archive. Historical
union-find roots are retained; observed-only family output avoids repeated scans
of all old families but does not bound checkpoint history. Resolve final-known
family roots with `canonicalize_stream_families.py`, retaining original IDs.

## Reproduce a rain-only reference window

`START_SNAPSHOT_INDEX` selects the preceding cumulative snapshot. N intervals
need N+1 consecutive files. First increments cannot be inferred without a
predecessor; accumulation resets, missing hours and grid changes fail.

```bash
python -u scripts/validate_real_data.py DATA_DIRECTORY \
  --wrf-cumulative --rain-kind cumulative \
  --wrf-time-source filename --wrf-grid-source attributes \
  --start START_SNAPSHOT_INDEX --hours 72 --crop 0 \
  --threshold 1 --bridge-radius 9 --min-size 1 --grid-km 4 \
  --workers 4 --id-parallel-mode frames --tau 0.35 --max-displacement 20 \
  --max-gap 1 --gap-tau 0.45 --gap-ambiguity 0.05 --event-overlap 0.1 \
  --gap-conflict-policy endpoint_overlap --gap-min-pixels 16 \
  --adjacent-policy overlap_ranked --velocity-reset off \
  --score-policy coherent_adjacent --family-map-scope observed \
  --chunk-frames 6 --stream --save-labels \
  --sequence-id UNIQUE_SEQUENCE --output-dir results/reference_run
```

Filename chronology does not independently verify the precipitation source year;
grid attributes do not verify geographic placement. Validate a coordinate reference
against this full-grid shape/projection before geographical analysis. Reference
coordinate year and rainfall year need not be equal. Do not silently mix climate
members or reset sequence identifiers when handing over checkpoints.

Branch IDs persist through continuation; events create new branches. Family roots
can change on later merges, so use namespace + branch ID for immutable branch keys,
and an additional final-known family root for retrospective family statistics.

Next scientific stage: retain this reference unchanged while testing identification
threshold/envelope sensitivity, rainfall conservation and original-STEP controls.
Evaluate untouched future and historical windows before committing the definition
to climatological production. If those tests require a logic change, make a new
version and rerun the fixed reference suite instead of silently rewriting v1.
