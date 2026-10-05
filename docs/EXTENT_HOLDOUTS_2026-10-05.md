# Cross-date extent verification, 2026-10-05

## Decision

Keep the named tracking reference frozen and its scientific defaults unchanged.
Keep 0.1-mm/h extent a separately named measurement experiment, not a certified
storm definition. The additional dates reproduce both its rain recovery and its
small-seed / large-extent tail. Neither greater coverage nor an extreme map is a
tracking-accuracy verdict. No minimum core area, higher core intensity or distance
cap was introduced.

Today's changes are validation, reporting and diagnostic plotting utilities;
`step/tracking.py`, `step/identification.py`, `step/envelope.py` and the reference
preset are not changed. Core reference commit remains `daf70111`; the extent
implementation under test was already saved at `94109b7`.

## Preregistered scope and stopping rule

Before examining window outcomes, select two full native-grid windows:

- July 15 00:00 cumulative predecessor; 72 hourly intervals ending July 15 01:00
  through July 18 00:00, 1996.
- August 15 00:00 cumulative predecessor; 72 hourly intervals ending August 15
  01:00 through August 18 00:00, 1996.

Both use the same 1749 x 2049 grid, 4-km spacing, hourly RAINNC differences,
filename chronology and projection metadata contract. No coordinates or internal
time arrays were inferred. These are held-out dates for extent-threshold
selection, but earlier seasonal engineering tests used the same year/model;
they are not independent observations, years, or future-climate validation.
June 1's existing 72-hour development run is shown alongside, not relabeled a
holdout. No date was selected because it yielded a desired result.

The requested cutoff was 14:00 Chicago (CDT, UTC-05:00), October 5. The bounded
runner records `PLAN.json` before launching any subprocess and reserves the last
10 minutes for reporting. It refuses existing output directories and stops its
own subprocess group on stage timeout/deadline. Child success contracts, not
exit status alone, determine completion. Partial output is retained for diagnosis
and is never marked successful.

## Coverage and connectivity

All extents reuse each date's exact saved A 1-mm/h seed rasters and tracking graph.
The paired B 0.5-mm/h identification/tracking variant is computed separately; it
is **not** the graph used by the frozen C extent comparison.

| Window start | Assigned rain: 0.5 extent | Assigned rain: 0.1 extent | Extra rain, percentage points | Multi-seed wet components: 0.5 / 0.1 | Max seed labels/component: 0.5 / 0.1 |
|---|---:|---:|---:|---:|---:|
| June 1, development | 87.733% | 95.552% | 7.819 | 126 / 558 | 4 / 16 |
| July 15 | 84.764% | 93.702% | 8.938 | 110 / 743 | 5 / 11 |
| August 15 | 88.278% | 96.067% | 7.788 | 150 / 792 | 3 / 20 |

The 0.1 extent also incorporates formerly unassigned >=0.5-mm/h rain connected
by weaker bridges: 261,747 pixels in July and 354,136 in August. These contribute
6.18% and 6.36% of added rain respectively. Thus lowering the extent threshold
is not exclusively adding a thin skirt below 0.5 mm/h. Multi-seed raw components
are grouping diagnostics; watershed still preserves separate seed ownership.
They do not certify false meteorological merges.

Previously assigned pixel ownership changes: zero in all three 72-hour windows.
This observation is distinguished from exact seed identity and nested support;
no general claim is made for other segmentation algorithms or seed choices.

## Seed-area tails and reviewed masks

| Window | 0.1 extent/core ratio p99 / max | 0.1 centroid drift p99, cells | Single-cell seeds' fraction of assigned rain | >=1000-cell seeds' fraction of assigned rain |
|---|---:|---:|---:|---:|
| June 1 | 325.67 / 11,848 | 28.90 | 0.146% | 84.349% |
| July 15 | 294.00 / 11,494 | 23.88 | 0.156% | 78.893% |
| August 15 | 463.97 / 8,248 | 28.45 | 0.110% | 88.003% |

Bins are descriptive existing seed snapshots, not a new scientific area filter.
Single-cell seeds can dominate extreme object ratios without dominating rain
mass. The 2-15-cell bin adds 1.236% of assigned July rain and 0.974% in August.
Large seeds still have broad extents: >=1000-cell seeds' p99 ratio is 8.11 in
July and 11.86 in August. Removing tiny seeds alone would not establish a
coherent physical storm definition, and its rainfall loss would need reporting.

Two selected maximum-area-ratio cases were replayed and visually reviewed with
unique target coloring, not repeating branch colors. Full-grid segmentation is
performed **before** cropping the figure:

| Interval end | Node / branch / local seed | Seed cells | 0.1 extent cells | Seed maximum mm/h | Extent mean mm/h | Drift, cells |
|---|---|---:|---:|---:|---:|---:|
| July 16 17:00 | 9521 / 6741 / 263 | 1 | 11,494 | 1.03345 | 0.21972 | 53.80 |
| August 15 04:00 | 991 / 783 / 191 | 1 | 8,248 | 1.00873 | 0.23169 | 68.77 |

The plots show broad weak-rain ownership attached to a seed only just above
1 mm/h; this is not evidence of a sustained convective core. The touched raw
0.1 components contain six distinct seed labels in the July case and four in
August; other seeds are still separately owned. Extent centroids do not feed
tracking. Reported drift is core-to-extent geometry, not motion; nominal cell
distance is not projection-corrected geographic distance.

## Integrity and execution evidence

Both holdouts pass all 72 seed hashes, catalog label/area checks, unchanged saved
graph sources, barriers, nested support and eligible-rain conservation. The
0.5 accounting/object geometry reproduces the earlier C calculation within each
paired run; unavailable old C raster hashes are not claimed verified.

- July A graph: 18,079 nodes, 5,271 edges; 17 split and 13 merge events, 9 gaps.
- August A graph: 16,218 nodes, 4,857 edges; 23 split, 17 merge and one complex
  event, 6 gaps.

An additional July execution replay uses 1 worker / 5-hour chunks versus
8 workers / 6-hour chunks. All eight comparisons passed: objects/IDs, full edges
and scores, events, final roots, timestamps, grid/units/cadence, identified raster
hashes and tracked raster hashes. This is an execution-invariance test, not a
speed benchmark. Results are in `july15/replay_w1_c5/summary.json` under
`a_reference_reproduction`; no branch renumbering was needed.
Plots and tests overlapped part of the runs, so no controlled scaling comparison
or RCC throughput projection is inferred from today's wall times.

Final full-suite verification: **160 passed, one skipped**. Existing NumPy /
netCDF4 binary-compatibility and xarray shape-deprecation warnings remain; no
test failed. New checks cover timezone-aware deadlines, no output overwrite,
timeout/failure completion guards, seed-area accounting and full-grid diagnostic
replay before cropping. Compilation checks passed for all three new utilities.

Validation and result analysis completed at **12:31 CDT**, before the requested
14:00 cutoff. Final repository recording/push follows this completed test run;
no additional scientific experiment is scheduled by this record.

## Reproduce and locate outputs

Run the fixed-window bounded suite from the experiment checkout (deadline must
be in the future and include an explicit timezone offset):

```bash
python -u scripts/validate_extent_holdouts.py /PATH/TO/1996_RAINNC \
  --output-dir results/extent_holdouts_new \
  --finish-before YYYY-MM-DDT14:00:00-05:00 --workers 8
```

The starts are deliberately fixed in this validation script; it is not a generic
date-selection optimizer. To test other preregistered dates, use the individual
commands in [the extent protocol](FROZEN_EXTENT_SENSITIVITY.md).

Local output root: `results/extent_holdouts_20261005_run1/`. Each window has
`prepared/`, `paired/`, `extents/`, per-stage logs and a selected `case_area/` map.
`comparison/REPORT.md` and `comparison/summary.json` include matched-date tables,
all seed-area bins and source artifact checksums. `PLAN.json`, `stage_status.json`
and `SUCCESS` retain the execution record. Data, rate cubes, catalogs and images
remain ignored by Git; only reusable utilities, tests and this record are shared.

Next scientific decision: whether the target is a general rainfall event or a
strong-core system, and how weak-only onset/decay should be represented. Today's
validation does not require that decision, but final identification certification
does. RAINNC alone does not establish convection, CAPE or updraft intensity.
