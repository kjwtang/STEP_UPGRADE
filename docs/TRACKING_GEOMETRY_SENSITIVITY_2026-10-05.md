# Rain definition and tracking sensitivity, 2026-10-05

## Scope and decision

The user's scientific focus is the effect of rain thresholds on tracking, not
changing catalog admission. Keep the complete existing catalog and the named
tracking reference unchanged. Strong-rain patch/phase diagnostics remain optional
sidecars; they neither remove objects nor feed matching, motion or event decisions.

Do not directly substitute a 0.1-mm/h expanded footprint into the reference
tracker. Its impact on lineage is substantial, not just a larger exported outline.
More associations, fewer short branches or more recovered rain do not establish
tracking accuracy. A 0.5 footprint also changes lineage and is not automatically
safe. These are isolated sensitivity runs, not a new scientific reference.

## Distinguish four experiments

| Variant | Identification | Tracking features | Exported extent | Interpretation |
|---|---|---|---|---|
| A, reference | 1 mm/h | Original 1-mm/h objects | Original | Frozen reference |
| B, threshold change | 0.5 mm/h | New 0.5-mm/h objects | Original B | Detector and tracker inputs both change |
| C, measurement only | Saved A | Saved A graph unchanged | A-seeded 0.5/0.1 partition | No tracking effect by construction |
| E, present experiment | Exact saved A objects/identities | A-seeded 0.5/0.1 partitions | Experimental graph features | Fixed entry, joint feature sensitivity |

Previous strong-phase tests and C extent tests cannot answer E's question. In
contrast, B changes the input objects: July's proximity-based 0.5-mm/h grouping
contains up to 81 distinct original A labels in one B group. This is not the count
of seeds in a raw connected weak-rain component. B also adds coreless weak objects,
so its total object count alone does not describe merging or fragmentation.

## Controlled E protocol

Reuse the existing June 1 development window and the preselected July 15 / August
15 windows: 72 hourly intervals each, 1996, native 1749 x 2049 grid, 4-km spacing.
Each rate input was prepared from cumulative RAINNC with filename chronology and
one predecessor snapshot. These dates are from one model/year, not independent
observational or future-climate validation.

`compare_tracking_geometry.py` writes its plan before processing. It runs three
independent state streams using the same frozen tracking parameters:

- `seed_control`: original 1-mm/h masks;
- `extent_0p5`: full-grid seeded watershed at 0.5 mm/h;
- `extent_0p1`: full-grid seeded watershed at 0.1 mm/h.

All seed pixels and frame/local-label identities must be retained, with no new
object entry. Coreless weak components remain unassigned; no distance cap or
stronger-core gate is introduced. Original source files and checksums are checked
before and after. Checkpoint/reload occurs every six hours. Deadline-aware runs
reserve ten minutes for recording before the user's 14:00 CDT cutoff.

The unexpanded control must exactly reproduce objects, complete edges/scores,
events, final family roots, chronology/grid/units and seed/tracked raster hashes.
Every experimental graph must pass forward-time, endpoint, event consistency,
monotonic nonrecycled ID and connected-branch checks. Graph differences use
immutable `(parent frame, local label, child frame, local label)` keys, not raw
branch numbers or repeating colors.

E changes overlap, rain-weighted centroid, area and mean intensity together.
Those fields in E's diagnostic `objects.csv` are expanded feature measurements,
not replacements for the authoritative seed catalog. This is not a causal
ablation of overlap alone. The frozen coherent adjacent score is
`max(0.60*IoU + 0.25*exp(-distance/radius) + 0.15*mean-intensity ratio)` over
complete raw and advected hypotheses. Changed geometry can change candidate
admission, score, competing assignment, events, velocity history and later gaps.
No scoring weights, event rules or thresholds were retuned for this experiment.

## Results: 216 distinct hourly intervals

Every original control passed all eight exact reproduction checks. Every E run
retained the same 51,066 seed-node snapshots across the three windows; all source
graph/input signatures remained unchanged. These are 216 distinct input hours,
not 648 independent hours merely because three geometries were processed.

| Window | Tracking geometry | Edges | Gap edges | Split / merge / complex events | Original seed rain associated with an incoming edge |
|---|---|---:|---:|---|---:|
| June 1 | Seed reference | 5,526 | 5 | 25 / 23 / 1 | 88.087% |
| June 1 | 0.5 extent | 6,860 | 18 | 59 / 48 / 6 | 89.312% |
| June 1 | 0.1 extent | 8,551 | 68 | 172 / 128 / 25 | 92.620% |
| July 15 | Seed reference | 5,271 | 9 | 17 / 13 / 0 | 84.900% |
| July 15 | 0.5 extent | 6,707 | 34 | 50 / 41 / 0 | 87.419% |
| July 15 | 0.1 extent | 8,693 | 97 | 206 / 193 / 21 | 91.097% |
| August 15 | Seed reference | 4,857 | 6 | 23 / 17 / 1 | 86.704% |
| August 15 | 0.5 extent | 6,196 | 29 | 58 / 56 / 1 | 88.260% |
| August 15 | 0.1 extent | 8,007 | 74 | 194 / 192 / 43 | 93.326% |

The association denominator excludes each window's first frame and always uses
the original seed catalog's precipitation sum. It does not include newly added
weak-rain mass and is not a correct-link probability. Across these windows, event
counts are 120 / 319 / 1,174 and gap-edge counts are 20 / 81 / 239 for seed / 0.5 /
0.1 geometry. These large shifts establish sensitivity, not false-event rates.

| Window | Extent | Added endpoint pairs | Removed | Retained pairs with event type changed |
|---|---|---:|---:|---:|
| June 1 | 0.5 | 1,558 | 224 | 49 |
| June 1 | 0.1 | 3,168 | 143 | 212 |
| July 15 | 0.5 | 1,607 | 171 | 43 |
| July 15 | 0.1 | 3,570 | 148 | 263 |
| August 15 | 0.5 | 1,491 | 152 | 55 |
| August 15 | 0.1 | 3,288 | 138 | 243 |

This is not exclusively a tiny-seed issue. Event counts whose **every original
endpoint seed** has >=1000 cells are June 30 / 42 / 55, July 17 / 27 / 54 and
August 28 / 51 / 87. At 4 km, 1000 cells is nominally 16,000 km2, but these
proximity-grouped seed objects need not be one contiguous convective core.

All variants' median branch span remains one hour. Single-snapshot branch
fractions decrease under E, but maximum branch span is not monotonically improved:
June is 56 / 56 / 51 hours and July 46 / 59 / 51 hours. A new event can legitimately
end a branch while its family persists. Neither more events nor a changed branch
number is, by itself, proof of an ID-recycling bug or physical discontinuity.

## Reporting and review contract

`summarize_tracking_geometry.py` uses identical original seed rain weights in
every graph's incoming-association fraction. Newly attached weak rain cannot
inflate that comparison's denominator or numerator. This remains an association
statistic, not correct-link precision/recall. Branch spans include gaps and event
cuts; 72-hour boundaries censor lifetimes.

Event and changed-edge endpoint areas are binned using original seed areas, not
expanded areas. This separates tiny-fragment sensitivity from changes also
affecting large seeds. These bins are descriptions, never new admission rules.

`plot_tracking_geometry_case.py` selects an added/removed pair by largest minimum
original endpoint rainfall, then an immutable-key tie break. It replays the
full-grid masks before plotting a crop. Selected extremes are not a representative
sample or an independent false-link verdict. Orange ownership and explicit seed
labels are used rather than cross-run branch colors.

Four full-grid cases were plotted and visually reviewed: July added pair
`(24,2)->(25,45)` and removed `(39,50)->(40,52)`; August added
`(39,62)->(40,80)` and removed `(37,84)->(38,62)`. The displayed large rain
regions change area and mean intensity substantially under E. For example,
August's added pair goes from 26,434 / 46,731 seed cells to 62,084 / 99,719 extent
cells at 0.1. Figures are consistent with geometry sensitivity but do not resolve
whether the added merge is correct. Saved incident edges expose competing links;
removing one endpoint pair need not entirely break the object's trajectory.

Incident-edge review distinguishes these cases. July's removed pair
`9037->9310` has no incident replacement in either expanded graph. August's
removed split pair `8696->8887` instead accompanies a retained parent link
`8696->8900` retyped from split to continuation. July's added split connects
`5175` to `5403` and `5446`; August's added merge connects `9099` and `9169` to
`9329`. These observations describe alternatives in the graph, not independent
physical judgments. A blanket "every removed pair is a broken storm" would be
incorrect.

## Verification and implementation boundary

Full regression suite: **172 passed, one skipped**. The six new regression tests
cover exact unexpanded control, expanded-node identity, alternate chunk/worker
execution on a synthetic sequence, source mismatch / corrupt seed hash,
deadline/no-overwrite guards, fixed seed rain weights, descriptive event bins,
deterministic case selection and incomplete plotting sources. All four real-data
case plots were also successfully replayed and visually checked. Existing
NumPy/netCDF4 compatibility and xarray shape-deprecation warnings remain.
Compilation and whitespace checks passed.

A full native-grid July replay changes execution from 4 workers / 6-hour chunks
to 1 worker / 5-hour chunks. All eight exact comparisons passed for **each** of
seed, 0.5 and 0.1 geometries: objects/IDs, edges/scores, events, final roots,
timestamps, grid/units/cadence and identified/tracked raster hashes. Differences
between the three geometries are therefore not explained by this execution
change. This replay repeats the same 72 hours; it is not another holdout or a
controlled speed benchmark. The check is saved as
`july_replay_w1_c5/EXECUTION_EQUIVALENCE.json`.

Scientific runs, accounting and exact execution comparisons completed at
**13:23 CDT**, before the requested 14:00 cutoff. Final regression and repository
recording follow; no further scientific parameter experiment is scheduled here.

Experiment code SHA256 is
`75b5ab7564ae04a6e1534080d943025c4728be06bb3955b08c911d484a59b186` in all three
registered plans. Local runtime: Python 3.13.5, NumPy 2.5.3, SciPy 1.18.1,
scikit-image 0.26.0, Matplotlib 3.11.2, xarray 2026.7.0, netCDF4 1.7.4.
The frozen core reference remains `daf70111`; no edits were made to
`step/tracking.py`, `step/identification.py`, `step/envelope.py` or the reference
configuration. Added utilities are diagnostic, not a new production workflow.

## Next bounded algorithm test, not implemented here

If an expanded overlap is worth exploring, test it as a secondary continuity
cue while retaining core-based motion/intensity and explicit competition/event
guards. Keep a pure reference control and a full-envelope comparator. Do not
silently replace all feature geometry, lower `tau`, or use stronger-core labels
to delete catalog nodes in the same experiment.

Review added, removed and retyped links, particularly new gap jumps and many-to-
many components. An expanded mask can stabilize a genuinely evolving rain region
but can also make distinct seed owners overlap across time. Neither explanation
follows from an edge count alone. Check a preregistered additional year/window
and cumulative-reset/month boundaries before any new reference promotion.

Frame-level stateless identification/geometry preparation may run in parallel;
tracking state mutation remains ordered. Independent years/members can run as
separate sequences. This sensitivity script is not a worker-scaling benchmark,
and it does not independently track spatial tiles/months and glue branch IDs.

## Reproduce

Use the installed experiment environment with the envelope extra. Run one
successful fixed-seed extent source at a time; choose a new output directory and
a timezone-qualified deadline in the future:

```bash
python -u scripts/compare_tracking_geometry.py PATH/TO/SUCCESSFUL_EXTENTS \
  --output-dir results/geometry_new --workers 4 --chunk-frames 6 \
  --finish-before YYYY-MM-DDT14:00:00-05:00 --timeout-seconds 1200
python scripts/summarize_tracking_geometry.py results/geometry_new \
  --output-dir results/geometry_report_new
python scripts/plot_tracking_geometry_case.py results/geometry_new \
  --output-dir results/geometry_case_new --geometry extent_0p1 --kind added
```

The local result root is `results/tracking_geometry_20261005/`, with `june01`,
`july15` and `august15` experiment directories. Results, raw data, arrays,
checkpoint archives and plots stay ignored; only small source/tests/docs belong
in the public-facing GitHub branch. The formal tracking implementation and preset
remain unchanged.
