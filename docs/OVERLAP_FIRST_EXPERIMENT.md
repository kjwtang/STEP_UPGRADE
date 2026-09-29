# Adjacent overlap-first experiment (opt-in)

Default `adjacent_policy='score'` is unchanged. This experiment is NOT a
scientifically validated release and does not force the seven RCC target pairs
to be linked. It is designed to test whether morphology-induced centroid jumps
should be allowed to break otherwise strongly overlapping adjacent objects.

## Rule

For candidates admitted by the existing candidate generator, require RAW
(unadvected) intersection >=16 cells, coverage of BOTH objects >=15%, and
coverage of at least one object >=50%. Recover the common raw intersection
from raw IoU and the two areas; never combine raw and predicted masks to qualify.
These constants are provisional experimental definitions, not literature-
calibrated thresholds or a claim that 16 cells has a universal physical scale.

Union qualifying pairs with existing event candidates. Multi-object components
produce split/merge/complex events first. Remaining isolated strong pairs are
continued before ordinary score-based assignment; other candidates keep normal
scoring. Gap scoring is unchanged. Actual original scores are retained in edge
outputs, so an accepted edge may be below tau. Weak baseline links are not
otherwise disabled; the event-union design can expand families and must be
reviewed for over-linking.

For adjacent continuations only, reset inferred velocity when area changes by
more than a factor of two or centroid displacement exceeds max_displacement.
The branch ID is retained. A reset is zero estimated velocity, not a measured
wind estimate. This is part of the opt-in experiment and is not separately
calibrated. Split/merge children already begin with no inherited velocity.

`adjacent_overlap.json` in each streaming chunk records accepted strong pairs
(including whether the score was below tau) and velocity resets. Candidate
scores are not probabilities. Checkpoints record the policy revision and reject
cross-policy resume. Ordinary baseline checkpoints remain compatible with the
default. The policy is supported by saved-stream replay diagnostics.

## RCC controlled run

Keep the same 72-hour input and gap guard as the previous experiment; add only
`--adjacent-policy overlap_first` and a new output directory. Do not overwrite
previous results or start from their final checkpoint.

```bash
/usr/bin/time -v python -u scripts/validate_real_data.py \
  /path/to/hourly/2005.07 --wrf-cumulative --rain-kind cumulative \
  --hours 72 --crop 0 --threshold 1 --bridge-radius 9 \
  --workers 4 --tau 0.35 --max-displacement 20 \
  --gap-conflict-policy endpoint_overlap --adjacent-policy overlap_first \
  --chunk-frames 6 --stream --save-labels \
  --sequence-id cimp_2005_july_full72 \
  --output-dir results/full_domain_72h_overlapfirst_run1

python scripts/compare_stream_tracks.py \
  results/full_domain_72h_gapguard_run1 \
  results/full_domain_72h_overlapfirst_run1 \
  --output-dir results/overlapfirst_comparison_run1
```

Comparison checks identical identification masks and object measurements, then
compares edges using (frame,local_label) endpoint keys, not branch numbers.
Review seven focus pairs AND all added/deleted/event-changed edges. More edges
is not success by itself. Look for large false event families and unnecessary
branch changes. Runtime/memory and full-domain chunk equivalence still require
real-data verification; synthetic checkpoint tests are not a substitute.

## Local regression coverage

- Persistent core with a transient distant lobe: adjacent ID continuity and
  velocity reset, rather than alternating gap chains.
- Split followed by merge: explicit event edges and fresh branch IDs.
- No overlap, tiny fragments and small-contact slivers: no overlap rescue.
- Truly empty middle frame: gap still permitted.
- Checkpoint restart: exact rasters, object/edge outputs and diagnostics;
  policy mismatch rejected.
- Default explicitly equals legacy score mode; existing baseline tests pass.
- Saved-statistics replay supports both policies with exact raster/CSV checks.

RCC 72-hour results recovered seven focus pairs but lost 30 other continuation
edges. Scientific validation remains pending; do not promote to default.

## Velocity-reset ablation

`--velocity-reset auto` (default) preserves historical behavior: morphology
resets for overlap_first, no resets for score policy. `off` disables all such
resets while retaining overlap-first rules; `morphology` explicitly enables
them. The effective policy is checkpoint-protected. Old overlap-first checkpoints
still replay with auto/morphology; they cannot resume with off.

Repeat the same run with `--velocity-reset off`, saving to
`results/full_domain_72h_overlapfirst_noreset_run1`. Compare with:

```bash
python scripts/compare_stream_tracks.py \
  results/full_domain_72h_overlapfirst_run1 \
  results/full_domain_72h_overlapfirst_noreset_run1 \
  --reference results/full_domain_72h_gapguard_run1 \
  --output-dir results/overlapfirst_noreset_comparison_run1
```

The reference section counts the original continuation pairs absent before,
lists which are recovered and their new event types, and lists any newly
missing reference pairs. Identification and object measurements must match all
three runs. Recovery is evidence about algorithm behavior, not ground truth.

## Paired diagnostic replay after ablation

```bash
python -u scripts/diagnose_velocity_ablation.py \
  results/full_domain_72h_overlapfirst_run1 \
  results/full_domain_72h_overlapfirst_noreset_run1 \
  --reference results/full_domain_72h_gapguard_run1 \
  --output-dir results/velocity_paired_diagnosis_run1
```

No NetCDF reads or identification reruns. Both full streams are replayed from
saved object statistics and masks; this is not independent verification of
rain ingestion or object measurement. Input masks/measurements must agree in
all three runs, and the two replay policies may differ only in effective
velocity reset. IDs are resolved independently from frame/local-label pairs.

Select all reference continue edges newly absent in reset-off, up to three
recovered controls, and all adjacent edges newly classified as events. In the
current 72-hour case this selects 10 lost pairs, 3 controls, and 2 split edges.
Report area, raw/advected intersection and coverage, velocity/history, score,
candidate admission, event eligibility, and actual incident edges. Raw overlap
is input evidence; advected overlap depends on prediction. A policy change can
alter earlier matches, so velocity differences are not necessarily attributable
to a reset at only the current frame.

Only after exact raster and all four catalog CSV checks succeed for both runs
are the top-level REPORT.md, paired_diagnosis.json and SUCCESS written. If the
second replay fails, partial subdirectories may remain; they are not a complete
paired diagnosis. Output directories are never overwritten. No algorithm
parameters are changed by this script and no automatic correctness verdict is
assigned to gained/lost edges or split events.
