# Core/envelope experiments — not a validated production release

Frozen production baseline: commit `7159898` (tracking revision 3).
Event-audit tooling: commit `1a8ec83`. Experimental code lives on
`codex/core-envelope-experiments`; main and production defaults are unchanged.

## Scientific interpretation

Detecting seeds separately from their extent is established methodology, e.g.
[tobac segmentation](https://gmd.copernicus.org/articles/12/4551/2019/) and
[PyFLEXTRKR detect-and-spread](https://gmd.copernicus.org/articles/16/2753/2023/).
The latter's cold-cloud thresholds are NOT precipitation-rate thresholds.
Neither reference validates our particular 1.0/0.1 mm/hr pair.

Here 1.0 mm/hr defines an operational rain seed, not a convective core. Requiring
such a seed excludes weak-only rain and may truncate developing/decaying storms.
Keep seed presence and rain extent separate from any eventual convective-storm
classification using stronger rain, updraft and environmental evidence. CAPE
alone is not proof of active convection. Do not impose a stronger threshold
without testing its dependence on grid spacing and hourly averaging.

## Two explicitly different interpretations

* `partition`: watershed on negative rain rate, seeded by all core pixels.
  Connected eligible weak rain is allocated to seeds, never across dry/NaN cells.
  Distinct seeds remain distinct even when their envelopes meet. Watershed
  boundaries and plateaus are algorithmic assignments, not observed boundaries.
* `system`: connected components of rain >= envelope threshold, retained only
  if seeded. A weak bridge can merge multiple seeds into one envelope. This is
  a deliberate sensitivity comparator, not an endorsement of unrestricted merges.

Default experimental seeds reproduce existing STEP radius-based identification;
`--seed-mode connected` removes that morphological joining and uses 8-neighbour
connected cores. Do not conflate changes in seed grouping with envelope growth.
`--min-core-cells` is configurable, but a one-cell seed is only the baseline
control, not scientific qualification. Convert area requirements to physical area
using actual grid information before production use. No spatial smoothing, dry
gap filling, distance cap, temporal hysteresis or HDBSCAN is silently introduced.
Experiments identify envelopes only: storm tracking remains on the frozen cores.

## Tuesday RCC run (separate checkout; do not overwrite existing results)

On login node, from the current repository:

```bash
git fetch origin
git worktree add ../STEP_envelope_experiments origin/codex/core-envelope-experiments
cd ../STEP_envelope_experiments
python -m venv .venv
source .venv/bin/activate
pip install -e '.[envelope]'
```

On compute node, activate that experiment environment and run:

```bash
python -u scripts/compare_envelopes.py \
  ../STEP_UPGRADE/results/tracking_revision3_600 \
  --output-dir results/envelopes_baseline_seeds
```

Only saved unthresholded rain rates are needed, not the original 2GB files.
Six variants (partition/system x .1/.2/.5 mm/hr) write frame-streamed label arrays
and CSV statistics to a new directory. Default seeds must exactly reproduce
the saved baseline. Core-size and connected-seed experiments can use additional
new output directories. Runtime includes all identification/I/O, not tracking.

Before promotion: inspect weak bridges, coreless excluded rain, envelope/core area
ratio, retained rain fraction, domain-boundary truncation, and sensitivity to core
size, seed grouping and thresholds. Compare historical/future climates using the
same definition. If later tracking envelopes, explicitly check centroid drift,
split/merge inflation and month/chunk equality; current core tracking tests do not
validate envelope tracking. Test full domain memory/runtime on RCC before scaling.

Status: local synthetic validation only; RCC and scientific validation pending.
