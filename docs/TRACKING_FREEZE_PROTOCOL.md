# Tracking freeze candidate: balanced overlap v2

Status: experimental, not yet frozen. Default scientific policies are unchanged.
This protocol concerns hourly 4-km precipitation objects, not a validated MCS or
convective-storm classification. Keep identification at 1 mm/h, bridge radius 9,
and minimum size 1 while comparing tracking. A 0.1-mm/h envelope is a separate
scientific experiment and must not silently redefine the tracked core.

## Fixed candidate

Use `--tau 0.35 --max-displacement 20 --max-gap 1 --gap-tau 0.45`
and `--gap-ambiguity 0.05`, with these opt-in policies:

```text
--gap-conflict-policy endpoint_overlap
--adjacent-policy overlap_balanced
--velocity-reset off
--score-policy coherent_adjacent
--gap-min-pixels 16
--family-map-scope observed
```

Balanced v2 uses stationary-mask intersection of at least 16 cells. Isolated
overlap continuations require at least 30% coverage of both endpoints. Event
topology has a separate rule: at least 15% of both endpoints and 50% of either.
This preserves unequal multi-way split/merge evidence before one-to-one
assignment. Ordinary score-qualified links remain possible. These constants
are experiment definitions, not universal physical thresholds.

The rejected v1 experiment applied the symmetric 30% rule to events as well;
a four-way split counterexample exposed that error. Its outputs must not be
presented as v2 validation. Checkpoints encode `adjacent_overlap_balanced_v2`.
The 16-cell gap floor excludes tiny endpoint gap links, not observed objects.

## Evidence required before freezing

1. Synthetic split, merge, four-way split, tiny fragment, true gap, competing
   assignment, monotonic/non-reused branch IDs, and checkpoint tests pass.
2. Performance-only changes preserve full rasters, catalogs and graph edges.
3. Fixed-policy independent windows include June/July continuity, a July week,
   and an August case. Compare additions **and** removals; baseline is not truth.
4. Review large-object interruptions and event masks visually. More links or
   fewer gaps do not demonstrate scientific accuracy.
5. Worker modes and alternate checkpoint schedules preserve results exactly.
6. Record commit, source hashes, input time/grid checks, timing, and limitations.

Filename timestamps are authoritative only with explicit
`--wrf-time-source filename`. Rain-only files without coordinates require
`--wrf-grid-source attributes`; grid attributes and dimensions must stay equal.
This does not independently verify geography. Accumulation decreases and missing
hourly predecessors fail; do not clamp resets into artificial rainfall.

## Parallelism and identities

`--id-parallel-mode frames` parallelizes independent identification frames.
Use chunks at least as large as the worker count. `tiles` parallelizes only
halo-aware dilation, then resolves connectivity and IDs globally; it cannot
introduce tile-seam storm IDs. Measure speed rather than assuming a gain.
Tracking state updates are causal and remain sequential. Do not concatenate
independently tracked temporal chunks and call that equivalent tracking.

Vectorized relabel/raster lookup removes repeated whole-domain scans without
changing first-pixel label order or branch allocation. Branch IDs are monotonic
within a sequence; a split/merge intentionally creates new branches. Family
roots can change after a later merge: that is distinct from renumbering a branch.
Observed family-map output avoids scanning all closed historical families every
frame, while preserving the complete union-find mapping in checkpoints. For
whole-sequence canonical family IDs, resolve every branch against the final
checkpoint; chunk catalogs represent roots known at that chunk's end.

## Scientific references

[PyFLEXTRKR](https://gmd.copernicus.org/articles/16/2753/2023/) motivates explicit
overlap lineage and configurable motion handling; its example thresholds and
cadence are not calibration for this dataset.
[tobac v1.5](https://gmd.copernicus.org/articles/17/5309/2024/) also distinguishes
feature tracking and split/merge handling. Neither reference makes STEP's
original raster IDs a ground-truth lineage catalog.
