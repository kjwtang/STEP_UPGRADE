# Change log

Commit hashes identify reproducible experiments; package version 0.3.0 alone
does not distinguish experimental revisions. Record `git rev-parse HEAD` with
each validation run. Do not update a checkout while its run is active.

## Unreleased — balanced v2 and scaling candidate

- Separate balanced continuation evidence from unequal split/merge topology;
  add explicit gap endpoint size floor. Scientific defaults remain unchanged.
- Vectorize stable relabel and branch raster construction; add exact halo
  dilation mode with global connectivity and observed-only family-map output.
- Stream rain-only hourly filename inputs with explicit grid-attribute checks;
  record source hashes. Add synthetic topology, seam and checkpoint tests.
- See `docs/TRACKING_FREEZE_PROTOCOL.md`; real-data v2 validation is pending.

## Unreleased — coherent adjacent-score experiment

- Opt-in `--score-policy coherent_adjacent`: compare complete stationary and
  constant-velocity scores instead of mixing overlap from either hypothesis
  with predicted distance. Multi-frame gap scoring and defaults unchanged.
- Checkpoint guards and replay diagnostics support the policy. This is an
  uncalibrated experiment; maximum-over-hypotheses may increase false links.

## Unreleased — full-grid local validation and overlap reuse

- Reuse exact child pixel sets and translated parent sets within each candidate
  pass. No persistent cache or score/threshold change. The 1996 June 72-frame
  candidate replay retains identical rasters, catalogs, scores and checkpoints.
- Local single-run tracking time fell from 68.14 to 25.86 seconds; total time
  from 98.27 to 55.42 seconds on a 1749-by-2049 grid with two ID workers.
  These timings exclude RAINNC conversion and are not a controlled RCC benchmark.
- Stream comparison now requires explicit focus pairs; historical node IDs are
  no longer silently applied to unrelated datasets. Replay can verify an
  alternate checkpoint schedule with `--checkpoint-every` and `--no-focus`.
- Add bounded RAINNC-only preparation using hourly filename timestamps. Internal
  time and coordinate arrays are optional for grid-space validation. Inputs
  remain read-only; missing files, grid changes and accumulation decreases fail.

## d86d4cd — replay performance and observability

- Replace the diagnostic duplicate full-frame candidate pass with requested
  pair-only admission checks; production full-frame assignment is unchanged.
- Reuse overlap evidence, index catalog rows by time, and report per-stage
  elapsed time with a 15-second heartbeat. Exact raster/catalog checks remain.
- Local tests verify pair-local candidate results against full-frame candidates.
  RCC speedup is not yet measured; interrupted outputs are not resumed.

## 413111d — paired velocity diagnostics

- Add exact paired replay of reset-on/off streams using frame/local-label
  identities. Automatically inspect newly lost reference continuations,
  recovered controls, and new event links. Export separate raw/advected
  coverage and motion history. Tracking behavior is unchanged.

## a57ae22 — velocity reset ablation

- Separate `--velocity-reset auto|off|morphology`; auto reproduces prior behavior
  and checkpoint configuration. Off disables resets without changing overlap
  thresholds. Effective-policy changes reject checkpoint resumption.
- Add optional reference-baseline recovery accounting to stream comparison.
- Test checkpoint/replay and identity continuity with reset disabled.

## e6efd0d — opt-in adjacent overlap-first experiment

- Add strong raw-overlap continuation/event policy, off by default, with
  morphology-triggered velocity resets and checkpoint configuration protection.
- Add saved-stream edge comparison by frame/local-label identity; support the
  policy in exact tracking replay. See `docs/OVERLAP_FIRST_EXPERIMENT.md`.
- RCC 72-hour comparison restored all seven focus pairs, but removed 30 other
  continuation edges without added same-endpoint replacements. Not promoted to
  default; velocity-reset ablation is needed before causal interpretation.

## 2a80a50 — adjacent pair diagnostics

- Add `replay_stream_pairs.py`: saved object-statistics/mask replay with exact
  raster and CSV checks before reporting adjacent-link rejection evidence.
  No algorithm change; source chunk outputs are read-only.

## cc85959 — repository hygiene and identifier regression

- Ignore local results, numerical inputs, logs and common credential files.
- Test persistent branch IDs despite local-label changes and earlier-object
  disappearance; test non-reuse after checkpoint resume.
- No tracking or identification behavior change in this maintenance update.

## 1641731 — opt-in gap conflict guard

- Experimental intermediate-object overlap veto, off by default; not a repair
  of adjacent links and not scientifically validated on RCC yet.
- Stream outputs include blocked-candidate diagnostics.

## 9be9522 — saved gap audit

- Inspect saved chunk labels and intermediate graph links; no algorithm change.

## c5f5b76 — progress reporting

- Cumulative tracked-frame percentage and final completion phases.

## 99240c0 / 657f2ca / 4532608 — research experiments

- Original STEP scientific context, overlap-event controls, core/envelope
  experiments. These do not establish a validated production configuration.

## 1a8ec83 / 7159898 — baseline tracking diagnostics and fixes

- Event audit, candidate gate changes and threshold-aware assignment.

## Public-release checklist

- Keep main as the baseline until reviewed experiments are ready to merge.
- Keep scientific inputs and generated results outside Git; publish only small,
  deliberately selected reproducible examples with documented permissions.
- Review tracked files AND Git history for credentials, personal paths, and
  identifying author metadata before making a private repository public.
- Resolve licensing and upstream attribution before distributing a release;
  this repository currently has no top-level LICENSE. Do not guess a license.
- Distinguish node, branch and family IDs. IDs are sequence-scoped; filtering
  reports must not renumber them. Split/merge can create new branches. Different
  reruns may allocate different IDs and require explicit correspondence.
- Consolidate historical experiment guides into a clear release guide after
  validation, without deleting evidence prematurely or rewriting history here.
