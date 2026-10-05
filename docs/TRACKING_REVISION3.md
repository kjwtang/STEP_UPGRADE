# Tracking revision 3: candidate gates and assignment fixes

Historical revision note. This change did not alter identification, rainfall thresholds,
score weights, tau=.35, gap thresholds, or event thresholds.

Candidates are the union of objects near the current centroid, near the predicted centroid,
and with substantive raw/advected mask coverage. Coverage uses `event_overlap`;
even a zero setting requires nonzero pixel intersection. A bounding-box-center KD-tree
conservatively screens possible overlaps, followed by actual pixel-coverage checks.
No all-pairs rain-pixel distance matrix is built. Bounding-box intersection alone cannot
admit a distant object. Scoring distance remains prediction error, isolating the candidate-policy change.

Subthreshold candidates are removed before assignment. Valid candidates are assigned by
connected component, with an unmatched option for every parent. The objective maximizes
total valid-edge score without requiring every object to match. Split/merge detection
still precedes one-to-one assignment; additional candidates may reveal previously missed
events, which require inspection of maps and edge tables.

`tracking_config` gained `algorithm_revision=3`. Earlier checkpoints are rejected;
restart the sequence with the new revision. Replay mismatch when newer code is applied to
old results is an intentional safeguard. Preserve the earlier diagnosis instead of overwriting it.

## RCC rerun without rerunning original STEP

Update the development checkout on a login node:

```bash
cd /project2/moyer/kjwtang/myproject/STEP_validation/STEP_UPGRADE
git pull --ff-only
```

Activate the separate upgraded environment on a compute node and run from the repository:

```bash
python -u scripts/recheck_tracking.py \
  results/original_vs_upgrade_600_mem24 \
  --output-dir results/tracking_revision3_600
```

Only upgraded tracking is rerun. Identification, hourly rainfall, and both original-tracking
groups are reused read-only through symlinks; retain their source directory. Tracking starts
from empty state. A separate six-frame chunk replay reloads JSON checkpoints and requires
exact agreement of all rasters and node/edge/event/family tables. New three-way comparison
plots and `edge_changes.json` include input hashes, source version, and parameters.
Compare immutable `(frame, local_label)` object identities, not branch numbers that may shift.

```bash
cat results/tracking_revision3_600/edge_changes.json
```

Inspect `focus_pairs`: does the object pair formerly labeled 36 -> 76 now have a continue
or event edge? The former 12 -> 26 pair remains below tau; this revision does not lower tau
to connect it. Review all added, removed, and changed-event edges, not only restored focus links.
Numerous new split/merge events require checking for false links between neighboring systems.

## Verified regression scenarios

- Motion reversal or centroid deformation breaks prediction, but the current-position route remains available.
- Both centroid routes miss, but substantive pixel overlap admits a candidate.
- Bounding boxes intersect without mask overlap, far from both centroids: no candidate is admitted.
- The 2x2 counterexample selects .90 rather than .80; all-ineligible scores return empty assignment.
- Split/merge, gaps, cross-chunk equivalence, and rejection of old checkpoints.

Real RCC effects and performance require the rerun above. Synthetic passes do not complete scientific validation.
