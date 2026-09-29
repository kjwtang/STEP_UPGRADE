# Change log

Commit hashes identify reproducible experiments; package version 0.3.0 alone
does not distinguish experimental revisions. Record `git rev-parse HEAD` with
each validation run. Do not update a checkout while its run is active.

## Unreleased — velocity reset ablation

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
