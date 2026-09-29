# Change log

Commit hashes identify reproducible experiments; package version 0.3.0 alone
does not distinguish experimental revisions. Record `git rev-parse HEAD` with
each validation run. Do not update a checkout while its run is active.

## Unreleased — repository hygiene and identifier regression

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
