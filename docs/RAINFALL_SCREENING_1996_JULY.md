# Frozen-reference rainfall screening: 1996 July

This is a diagnostic of the **1-mm/h / radius-9 reference on the full grid**,
not observational storm truth and not a reason to silently change the definition.
It uses the saved 744-hour July control catalog: 175,329 object snapshots,
53,732 accepted edges and 122,383 observed branches. Family roots are resolved
to the last checkpoint before retrospective grouping.

## Counts and labeled-rain contribution are different diagnostics

| Object area (grid cells) | Object snapshots | Share of snapshots | Share of labeled rain |
|---|---:|---:|---:|
| 1–15 | 93,145 | 53.13% | 0.298% |
| 16–99 | 41,322 | 23.57% | 1.530% |
| 100–999 | 27,361 | 15.61% | 11.380% |
| >=1,000 | 13,501 | 7.70% | 86.791% |

“Labeled rain” is the sum of rain rate over identified >=1-mm/h grid pixels,
weighted equally by grid cell, accumulated across catalog snapshots. It does
not include subthreshold rain and is **not geographic area-corrected volume**.
Those fractions cannot be used as the fraction of all rainfall captured.

The full 2208-hour JJA control gives a similar separation: 269,496 of 504,350
object snapshots (53.43%) are smaller than 16 cells, representing 0.269% of
labeled rain. This is another diagnostic of the same 1996 dataset, not an
independent year or a historical/future validation sample.

About 84.88% of branches have only one observed snapshot. The median observed
branch span is one hour and the maximum is 142 hours, gaps included. Split/merge
events intentionally terminate/create branches; a short branch is not necessarily
a broken meteorological system. Long rain-band tracks are not independently
certified convective storms either. Time and domain boundaries censor tracks.

Excluding the first frame when computing incoming associations, 30.59% of object
snapshots have an accepted incoming edge, representing 85.67% of labeled rain.
For >=1,000-cell objects, the corresponding fractions are 85.83% and 88.88%.
These are **graph association coverage**, not accuracy: genuine genesis/lysis
may be unlinked, and accepted links may be wrong. There is no human/observational
reference-label set supporting a precision/recall claim.

## Implications for the next identification experiment

1. Keep the frozen reference unchanged while running a separately named sweep.
   Do not remove small objects merely because they dominate counts: some can be
   real isolated cells. Rain-volume relevance and physical classification are
   different objectives.
2. Audit moderate/large object's interruptions by masks, raw/advected coverage,
   independent timestamps, and birth/decay context, not all-node link ratio alone.
   Size bins and labeled-rain weighting help prioritize inspection; they do not
   supply ground truth.
3. Before using 0.1-mm/h envelopes, measure **all-rain** accounting directly from
   unthresholded hourly inputs: seeded assigned rain, coreless weak rain, missing
   cells and unassigned components. None of those quantities can be reconstructed
   from >=1-mm/h catalogs alone.
4. Separate an operational rain seed from a physical convective core. Requiring
   stronger/long-lived seeds, or discarding small branches, changes which
   developing/decaying and weak-only systems are represented. Rain-only fields
   do not establish updraft or CAPE evidence.
5. Keep the existing core-tracking graph fixed for an initial envelope study.
   A changed envelope centroid or weak bridge can alter subsequent tracking;
   testing envelopes alone does not validate tracking them.

Reproduce this screening on a successfully saved stream:

```bash
python scripts/summarize_reference_catalog.py results/REFERENCE_STREAM \
  --output-dir results/REFERENCE_SCREENING
```

The script saves size/lifetime plots and statistics locally, including the
saved identification threshold when available. It does not upload source rain,
rasters or generated catalogs to GitHub. Existing
[core/envelope definitions and caveats](CORE_ENVELOPE_EXPERIMENTS.md) remain
separate from the frozen tracking reference.

## Preliminary unthresholded-rain accounting: June 1, 72 hours

A separate full-grid accounting experiment uses the unthresholded hourly
RAINNC increments, not the July labeled catalog. Equal grid-cell rain-rate
sums across the 72 frames give:

| Threshold (mm/h) | Share of all provided valid RAINNC rain | Eligible share of valid grid-time cells |
|---:|---:|---:|
| 0.1 | 97.15% | 8.71% |
| 0.2 | 94.77% | 6.61% |
| 0.5 | 88.47% | 4.13% |
| 1.0 | 79.75% | 2.58% |
| 2.0 | 66.84% | 1.41% |
| 5.0 | 45.20% | 0.52% |

The 0.1-mm/h eligible footprint is approximately 3.38 times the 1-mm/h
footprint in this window. This establishes a potentially important weak-rain
contribution, not that every eligible pixel belongs to a seeded storm.
RAINNC-only accounting excludes any unavailable RAINC contribution and is
not map-area-corrected volume. Raising the seed threshold cannot be justified
from rain contribution alone either.

On 12 sampled frames (every sixth frame), existing partition/system envelopes
at 0.1 and 0.5 mm/h passed exact 1-mm/h reference core-mask equality and equal
assigned-support checks between envelope modes. This does not certify weak
bridges, ownership, convective cores or tracking of expanded masks. The
reference tracking and identification parameters remain unchanged. This is
an exploratory sample from the same 1996 year, not independent validation.

Local reproducibility evidence is saved under
`results/threshold_accounting_run1/summary.json`; raw inputs and experimental
arrays are intentionally excluded from GitHub. A next scientific experiment
must inspect coreless weak rain, multicore bridges, ownership, genesis/decay
and independently selected windows before promoting an envelope definition.
