# Tuesday experiment matrix (saved-data first)

No RCC runs have been performed for these additions. Main stays unchanged.

| Version | Git reference | Purpose |
|---|---|---|
| T3 | `7159898` | Baseline candidate/assignment fixes; already run on RCC |
| Audit | `1a8ec83` | All-frame tau sensitivity tool |
| E1 | `4532608`, `codex/core-envelope-experiments` | First partition/system seed-envelope experiment |
| E2/T4 | `codex/literature-informed-experiments` | Physical controls/statistics and optional event policy; unvalidated on RCC |

Use a separate worktree/environment on the login node, leaving the old installation
and output directories intact. From the existing STEP_UPGRADE repository:

```bash
git fetch origin
git worktree add ../STEP_literature_experiments origin/codex/literature-informed-experiments
cd ../STEP_literature_experiments
python -m venv .venv
source .venv/bin/activate
pip install -e '.[envelope]'
```

On an allocated compute node, activate that environment and run:

```bash
python -u scripts/run_tuesday_suite.py \
  ../STEP_UPGRADE/results/tracking_revision3_600 \
  --output-dir results/tuesday_suite
```

The runner stops on failure, saves stage logs and `suite_status.json`, and writes
top-level SUCCESS only when everything completes. Default timeout is 900 seconds
per stage; change explicitly for larger cases. It does not impose an RSS limit:
stay inside the allocated Slurm memory and use the 600-cell case first.

Runs:

1. **T3 control and tau sensitivity:** .35 (exact saved raster replay), .30, .325.
2. **Matched event control:** score-and-overlap, overlap .5, minimum 2 shared cells.
3. **T4:** overlap-only events with the same .5/2; continuation/gap rules unchanged.
4. **E2 baseline seeds:** partition/system each at .1/.2/.5, exact baseline seed check.
5. **E2 connected seeds:** same six envelope settings, without original morphology.

Items 1-3 share the `events` stage. Items 4-5 are separate envelope stages. Twelve
int32 envelope cubes plus two core cubes require about 484 MB for 24x600x600,
excluding CSVs and logs (about 461 MiB). Check disk capacity before full-domain use.
Frame streaming bounds working arrays, not the total output size.

The default suite intentionally does NOT invent grid spacing, minimum physical
core area or a growth radius. Once spacing is verified, add an optional bounded
stage with `--spacing-km DY DX --max-distance-km R` (replace tokens with verified
values and a documented test radius). This is nominal Cartesian distance, not
geodesic distance. It is distance to any retained core, not the assigned core.

Physical area experiments use `compare_envelopes.py` directly with
`--min-core-area-km2 A` plus either `--cell-area-km2 CELL_AREA` or
`--cell-area-file areas_same_crop.npy`. Provide actual physical areas or explicitly
label nominal-area approximations. Do not reuse a cell-area array from another crop.

## Outputs to review together

* `events/summary.json`, `events/edge_changes_vs_baseline.json`: edge/event counts
  and added/removed/changed-event edges keyed by frame/local-label, not branch ID.
* `events/*/candidates.csv`: why event candidates pass/fail; only adjacent active
  candidates, not a complete audit of dormant-gap rejection.
* `envelopes_*/core_labels.npy`: retained seed labels, saved independently.
* `envelopes_*/*/objects.csv`: per-object intensity, area ratio, censoring flags and
  optional area-weighted volume rates (blank physical values without explicit areas).
* `envelopes_*/*/core_membership.csv`: many-to-many seed/envelope associations.
* `envelopes_*/*/statistics.csv`: frame-level coverage/rate sums. Legacy column
  `unassigned_weak_cells` counts all eligible unassigned cells >= lower threshold,
  not just cells strictly below the core threshold; use that interpretation.

No minimum event count, longest-track target or automatic “best configuration”
selection. Spatial review and scientific definition come before promotion.
See [method review](STORM_METHOD_REVIEW_2026-09-26.md) for sources and limitations.
