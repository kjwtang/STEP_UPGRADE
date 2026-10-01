# Mutual-dominant overlap experiment

Opt-in `--adjacent-policy overlap_ranked`. Defaults and balanced-v2 results
are unchanged. This is a testable relational ranking rule, not HDBSCAN and not
a probability model or validated physical storm classifier.

Keep balanced-v2 event topology and symmetric strong continuations. Before
ordinary score assignment, admit an otherwise unused pair when:

- Stationary-mask intersection is at least 16 cells and covers at least 15%
  of both endpoints.
- The pair has each endpoint's largest raw intersection across **all** its
  candidate neighbors, and at least twice the runner-up intersection. Ties reject.
- Neither endpoint was consumed by a significant event or strong continuation.

Candidates consumed by earlier stages still participate in ranking. Do not
remove a competing best match and then promote a weaker residual as dominant.
Only actual ranked continuations are marked `mutual_dominant_raw_overlap` in
the audit log. Motion-only overlap cannot establish this rescue condition.

This tests whether ranking can preserve moderate-overlap continuity without
lowering the global score cutoff or rainfall threshold. The geometric floors
and dominance factor are explicit experimental definitions, not literature
calibration. They must be checked for false connections and changed downstream
motion history, not just increased edge counts.

Fixed validation: reuse June-1 development and June-30 diagnostic windows;
evaluate untouched July-10 and August-20 72-hour windows. Do not call the reused
windows independent holdouts. Compare against balanced-v2 with identical
identification, and verify seven-frame checkpoint restoration. Preserve every
version by commit, rather than replacing old result directories.
