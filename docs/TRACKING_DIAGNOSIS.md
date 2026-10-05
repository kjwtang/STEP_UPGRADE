# Broken-link diagnosis without algorithm changes

Historical guide: replay old results only with the code that generated them.
After upgrading to 0.3, first run `recheck_tracking.py` as described in
[revision 3](TRACKING_REVISION3.md); do not replay old results with newer code directly.
The branch IDs below refer to the original diagnostic sample and may differ in later runs.

Input is the completed `original_vs_upgrade_600_mem24` output directory.
The script reuses saved rainfall, new identification labels, upgraded and hybrid tracking
rasters, and parameters. It does not read the original 2 GB files or rerun original STEP tracking.

```bash
python -u scripts/diagnose_tracking.py \
  results/original_vs_upgrade_600_mem24 \
  --output-dir results/diagnosis_600
```

The script replays all 24 upgraded frames and requires exact frame-by-frame raster equality.
It stops on mismatch so that different code/parameters cannot explain historical results.
By default it audits every adjacent object pair at 01:00 and 05:00, highlighting branch
12 -> 26 at 00:00 -> 01:00 and 36 -> 76 at 04:00 -> 05:00.
Use `--frames 1 2 3 4 5` for other transitions; these are zero-based child-frame indices.

Outputs are `diagnosis.json`, `REPORT.md`, and `pair_decisions.csv`. Records include raw
centroid displacement, prediction error, gates, raw/advected IoU, coverage, intensity ratio,
score, candidate admission, accepted edges, and other edges at the endpoints.
Scores outside candidate gates are counterfactual: they do not imply the pair entered assignment.

Decision categories:

- `outside_all_candidate_gates`: no candidate entry route passed; formerly `outside_prediction_gate`.
- `below_score_threshold`: admitted but below tau.
- `endpoint_used_by_other_edge`: eligible, but an endpoint was assigned elsewhere; inspect that edge before calling this a bug.
- `eligible_but_not_assigned`: eligible but unassigned; inspect competition and assignment results.
- `accepted_*`: an actual edge, including continue/split/merge. A branch-ID change need not be a broken link.

Original same-ID relations are comparison evidence, not ground truth; shared IDs can imply many-to-many relations.
The script also includes a 2x2 synthetic counterexample with tau=.35 and scores [[.90,.34],[.80,0]].
Assignment before threshold filtering chooses .80 instead of feasible .90; revision 3 filters first
and is expected to choose .90. This proves a historical assignment issue, not its causal role
in the two real-data breaks. The diagnostic script does not change thresholds, movement gates,
or production tracking code.
