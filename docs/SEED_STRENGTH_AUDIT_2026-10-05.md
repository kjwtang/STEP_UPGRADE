# Strong-rain phase sensitivity on the frozen graph

## Scope and preregistration

This is an audit of existing 1-mm/h objects, **not** alternative 2/5-mm/h
identification or tracking, a new convection classifier, or a default change.
Use the same June 1 development and July 15 / August 15 evaluation windows as
[the fixed-seed extent verification](EXTENT_HOLDOUTS_2026-10-05.md). Each contains
72 hourly intervals on the full 1749 x 2049 grid. The three windows are different
dates of the same model/year, not independent observations or future-climate
validation. No window was added or removed to obtain a desired conclusion.

Before reading strength outcomes, fix these sensitivity choices:

- Raw hourly rain rate: >=1, >=2, >=5 mm/h.
- Largest contiguous 8-neighbor strong patch inside the saved 1-mm/h object:
  >=1, >=4, >=16 grid cells. No dilation is applied to strong-patch area.
- Retrospective branch-phase qualification: at least 1 or 3 consecutive observed
  hourly frames with the chosen rate/patch condition.

At nominal 4-km spacing, the patch counts correspond to 16, 64 and 256 km2 in
grid space, not projection-corrected physical area. These numbers span a
diagnostic sensitivity matrix; they are not externally certified storm criteria.
An hourly output meeting a condition is not proof that rain exceeded it at
every instant within that hour.

The literature supports considering intensity, area and life history separately;
it does not prescribe this matrix for our RAINNC-only data. See
[the definition review](STORM_DEFINITIONS_REVIEW_2026-10-02.md). A 1-mm/h seed,
or even a 5-mm/h rain patch, does not independently establish convection.

## Two distinct selections

`snapshot` keeps a saved object's measured extent only at the times its strong
patch condition is met. It makes no duration requirement.

`branch_phase` qualifies an existing branch if its longest consecutive qualifying
run reaches 1 or 3 hours, then includes **all observed object snapshots in that
branch**, including its weaker phases. This is retrospective classification:
future observations can determine qualification of an earlier snapshot. It does
not change online tracking IDs, propagate qualification through split/merge to
another branch, or label the entire connected family a strong storm.

Gaps interrupt consecutive qualification. Event-ended branches and audit-window
boundaries truncate the assessed duration. Weak-only onset/decay stages outside
the original 1-mm/h graph are still absent. A three-hour branch condition must
not be presented as a measured minimum physical storm lifetime.

Rain accounting sums each retained node's **already saved** 0.5 or 0.1 extent
rain. Excluded ownership is not redistributed to surviving seeds. Reidentifying
at 2/5 mm/h and rerunning watershed would change seeds and potentially topology;
today's retained fractions are not predictions of that different experiment.

## Integrity and provenance

The auditor replays each 1-mm/h seed frame, checks all hashes and catalog areas,
and measures raw strong support. Counts and largest contiguous patches must not
increase when the intensity threshold rises. The 1-mm/h / one-cell snapshot
control must reproduce all original assigned extent rain. Both extent catalogs
must contain exactly the saved graph's node IDs and conserve their rain totals.
Graph, input, metadata and extent source files are checked unchanged after audit.

June's initial attempt stopped at a manifest mismatch before processing frames.
Read-only inspection found all 86 recorded signatures unchanged, no removed
files, and only two later manifest-coverage additions: `frame_accounting.csv`
and `c_envelope_objects.csv`. Compatibility is restricted to those exact two
paths; changed/missing signed files or new graph files still fail closed. The
final audit signs the expanded 88-file manifest in its own `PLAN.json` and
explicitly records these coverage additions. Old unsigned coverage is not
retroactively claimed verified.

The final three runs use one identical auditor source hash, stored in their
plans. No tracker, identification implementation, envelope implementation or
reference preset is changed. Stage concurrency and tests overlap, so their
wall times are not a controlled worker benchmark.

## Results

All three final runs passed: unchanged seeds, graph and extent sources; the
one-cell / 1-mm/h control preserves 100% of assigned rain. Source SHA-256 is
identical across runs:
`6a0f0918e3a22adfa68fa6d95ffc6b363adc4d1a5dabe20049ce9dad386dc08b`.

The following fractions are percentages of **original assigned 0.1 extent rain**,
not of all supplied rain and not accuracy scores:

| Rate / patch cells | Selection | Consecutive hours | June 1 | July 15 | August 15 |
|---|---|---:|---:|---:|---:|
| 1 / 1 | branch phase | 3 | 81.365% | 84.710% | 86.296% |
| 2 / 1 | each snapshot | -- | 97.515% | 97.064% | 97.633% |
| 2 / 4 | branch phase | 1 | 97.693% | 97.212% | 97.721% |
| 2 / 4 | branch phase | 3 | 80.399% | 83.334% | 85.130% |
| 5 / 16 | branch phase | 1 | 89.765% | 86.293% | 92.440% |
| 5 / 16 | branch phase | 3 | 72.809% | 74.507% | 80.379% |

For the 2-mm/h / four-cell / one-hour branch-phase choice, retained object
snapshots are 7,865 of 16,769 (June), 8,010 of 18,079 (July), and 7,159 of 16,218
(August). Those are 44-47% of snapshots but 97.2-97.7% of previously assigned
0.1 rain. Relative to **all supplied rain**, coverage is instead 93.347%, 91.089%
and 93.877%. This suggests that an optional descriptive strong-phase tag can
separate many low-rain snapshots without replacing the complete rainfall catalog.
It does not validate this particular threshold or prove the rejected objects
unimportant for event counts, morphology, extremes or climate changes.

The three maximum-area-ratio single-cell cases reviewed in the extent work
(June node 7050, July node 9521, August node 991) each have just one observed
branch snapshot and zero >=2-mm/h pixels. A 2-mm/h phase tag excludes those
examples' owned extent from its selected subset. These are selected extremes,
not an unbiased false-object sample, and their exclusion is not threshold
certification.

The duration condition has a substantially larger effect than a modest rate
tag in these windows. Even the 1-mm/h / one-cell control loses 13.704-18.635%
of assigned rain when only branches with three consecutive observed hours are
retained. Do not quietly use that as a minimum physical storm lifetime.

### Duration-context flags

Each entry below is a percentage of original assigned 0.1 rain belonging to
branches rejected by the **1-mm/h / one-cell / three-hour branch control**.
Subsets overlap and must not be added:

| Rejected subset | June 1 | July 15 | August 15 |
|---|---:|---:|---:|
| All short consecutive observed branches | 18.635% | 15.290% | 13.704% |
| Touching audit-window boundary | 3.007% | 0.859% | 0.915% |
| Touching split/merge/complex endpoint | 8.942% | 0.664% | 3.493% |
| Containing a gap edge | 0.000% | 0.027% | 0.030% |
| Without those flags | 8.490% | 13.756% | 9.608% |

June is particularly event-sensitive: 47 rejected branches touching event
endpoints account for 8.942% of assigned rain. This records an association with
algorithmic segmentation, not a causal attribution or a proof that connecting
those branches as one lifetime would be correct. Unflagged rejection may still
reflect score/gate fragmentation or genuine short rain episodes; it is not
automatically physical short life. Window contact is possible censoring, not a
confirmed birth or decay.

### Decision

Keep the complete 1-mm/h rainfall lineage and preserve the original extent
budgets. Strong-phase diagnostics may be layered on it, but no selected cutoff
becomes an identification default or a convection/MCS label. A branch can contain
widely separated raw pieces because STEP's existing proximity grouping is still
in force; qualification by one strong patch does not establish the physical
coherence of all its weaker extent. Existing outside-graph weak stages are not
recovered by retrospective selection.

Before making strong-phase qualification a catalog-entry requirement, decide
whether the scientific target is general rainfall events or a strong-core subset,
then assess that definition on an untouched year/sequence. Real 2/5-mm/h
identification, weak-only temporal attachment and family-level qualification are
separate future experiments. No such change was made today.

### User decision and later data processing

After reviewing the audit, the user explicitly chose **not to change catalog-entry
conditions**. Keep the complete existing rainfall catalog. The sensitivity table
is a hypothetical accounting selection, not a filtered catalog or an implemented
classification rule. Further work may plan data handling without silently
activating any strong-core gate.

For later multi-year data, the proposed arrangement is:

1. Normalize cumulative RAINNC to verified hourly intervals, with explicit
   predecessor snapshots, missing-hour/reset handling and source provenance.
   Monthly/yearly file boundaries must not restart the accumulator implicitly.
2. Preserve the named 1-mm/h identification and ordered tracking-state stream,
   handing checkpoints across contiguous chunks/months. Independent members or
   climate sequences may run separately; within-sequence chunks are not
   independently tracked and glued by renumbering.
3. Retain baseline node/branch/event IDs and the complete rain accounting.
   Record alternative 0.5/0.1 extent budgets separately; expanded centroids do
   not feed tracking and omitted weak-only rain remains explicit.
4. Store numerical strong-rain diagnostics in a separate sidecar keyed by
   immutable sequence/frame/node identities: raw peak rate, threshold-pixel
   counts, largest contiguous patch, and observed branch-phase duration/context.
   Keep rates/areas available rather than baking a binary convection label into
   the catalog. Retrospective labels must be distinguished from online data.
5. Assess historical/future definitions identically, with full versus qualified
   rainfall budgets reported together if a later scientific classification is
   authorized. Do not infer future thresholds from whichever subset looks most
   favorable in these three 1996 windows.

This is a processing plan, not a newly implemented production pipeline. Current
audits already retain per-node numerical diagnostics; no sidecar is consumed as
a catalog filter and no automatic classification is enabled.

## Verification and timing

Final full test suite: **166 passed, one skipped**. Existing NumPy/netCDF4
binary-compatibility and xarray shape-deprecation warnings remain; no test failed.
New tests cover disconnected patches sharing an object ID, noncontiguous IDs,
cross-object patch rejection, gap/weak-hour duration breaks, retrospective phase
selection, overlapping context flags, narrow legacy-manifest compatibility,
rain control, output overwrite, deadline reserve and invalid-source guards.
The audit script compiles successfully. Final three runs and analysis completed
at **12:49 CDT**, before the continuing 14:00 limit. No background scientific
experiment is left running by this record.

## Reproduce

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python -u scripts/audit_seed_strength.py SUCCESSFUL_FIXED_EXTENTS \
  --output-dir results/seed_strength_new \
  --finish-before YYYY-MM-DDT14:00:00-05:00
```

Deadline is optional for reuse; when supplied it must include an explicit UTC
offset. The audit reserves the final ten minutes for reporting, never marks a
partial run successful, and refuses output overwrite. Source rate cubes and
reference catalogs must still be accessible at their recorded paths.

Outputs: `PLAN.json` (fixed choices, full signatures, auditor source hash),
`object_strength.csv` (per-node/threshold strong counts and raw patch properties),
`summary.json` (all 27 selections plus three-hour branch context), `REPORT.md`
and `SUCCESS`. Window/event/gap context flags overlap and must not be summed;
they are not causal labels of short physical storms.

Local final output root: `results/seed_strength_20261005_final/`, with `june01/`,
`july15/`, `august15/`. Initial diagnostic outputs remain separately under
`results/seed_strength_20261005/`; the June attempt did not create a successful
run there. Raw data and generated catalogs remain Git-ignored.
