# Paired rainfall-threshold sensitivity — experimental, defaults unchanged

Measured first-window findings are saved in
[1996 June 72-hour results](RAIN_THRESHOLD_JUNE72_RESULTS.md).

The frozen `rainfall_lineage_4km_1h_reference_v1` remains unchanged. This study
changes the rainfall definition, not the tracking scores, gates, gap guard,
velocity reset, minimum size or radius. A better-looking or longer track is
not independent evidence that the new definition is meteorologically correct.

## Three deliberately different variants

| Variant | Identified / tracked pixels | Reported rain extent |
|---|---|---|
| A | >=1.0 mm/h, radius-9 STEP grouping | Same pixels |
| B | >=0.5 mm/h, otherwise identical ID/tracking settings | Same pixels |
| C | Exactly A's seeds and tracking graph | Connected >=0.5-mm/h rain partitioned among A seeds |

Rates are hourly RAINNC increments on the same 4-km grid. Minimum object size
remains one cell in both A and B, deliberately isolating the threshold change.
The radius stays nine cells: it is a dilation radius, not a maximum separation.
Lowering the threshold can change morphological group connectivity as well as
rain extent. No threshold-specific tracking tuning is performed here.

C uses the existing marker watershed, with all A seed pixels as markers and
dry/NaN cells as barriers. Its envelope label maps back to the A branch ID for
that frame; there is no new C tracking call, score, event or lifecycle claim.
An A seed may have disconnected pixels under STEP grouping. The same branch
can therefore own disconnected C pieces, consistently with the existing seed
definition. A multi-seed weak connected component is saved as a bridge
diagnostic; it is not silently substituted for C partition ownership.

Separating feature detection from extent segmentation has methodological
precedent in [tobac](https://gmd.copernicus.org/articles/12/4551/2019/), but that
paper does not validate the 1.0/0.5-mm/h pair. Rain-only seeds are operational
objects, not independently diagnosed convective cores. Weak-only rain and
birth/decay periods can be excluded by a seed requirement; report that loss.

## Reproduce on explicitly prepared hourly rates

Do not pass raw cumulative RAINNC as the NPY cube. The rate preparation needs
one extra predecessor snapshot, a verified hourly filename cadence and explicit
reset handling. Source files stay read-only. First prepare a bounded window:

```bash
python -u scripts/prepare_rain_only.py RAINNC_DIRECTORY \
  --start 1996-06-01T00:00:00 --hours 72 \
  --output-dir results/june72_prepared
```

Then, on the same pinned checkout and with the envelope extra installed:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python -u scripts/compare_rain_thresholds.py \
  results/june72_prepared/rain_mm_h.npy \
  --metadata results/june72_prepared/input_metadata.json --units mm/h \
  --workers 8 --chunk-frames 6 \
  --output-dir results/threshold_paired_june72
```

Optional `--reference FROZEN_A_STREAM` requires exact reproduction of its
complete A catalog, scores, events, final family roots and every raster hash.
Use `--no-plots` to suppress four sample full-grid maps. The new experiment
never overwrites an output directory and does not implement partial-run resume.
Its `frame_hashes.json` files are diagnostic records, not transactional
`run_reference_stream` commit markers. Root `SUCCESS` means all paired checks
and report generation completed, not that B or C should be promoted.

## Comparison outputs and interpretation

- `A_1mm/` and `B_0p5mm/`: complete chunk catalogs, checkpoints, frame hashes
  and structural audits. Full raster archives are not written by this study.
- `frame_accounting.csv`: all supplied rain, labeled A/B rain, seeded C rain,
  lower-threshold rain unassigned to C, cell counts and grouping counts.
- `a_to_b_memberships.csv`: pixel intersections of A seed labels with B labels.
  Numeric labels/node IDs cannot be equated across different definitions.
- `projected_edge_changes.json`: A/B edge comparison only where both endpoints
  have bijective shared-core object correspondences. Multicore absorption and
  weak-only endpoints are excluded and counted. These are not false-link rates.
- `weak_bridge_components.csv`: weak connected components containing multiple
  A seeds, preserving their membership relation. Additional radius-based B
  joining is reported separately from connected weak components.
- `c_envelope_objects.csv`: core/extent size ratio, retained A node/branch ID,
  rain sum and intensity-weighted centroid drift. Drift measures a possible
  consequence of future envelope tracking, not a change to C's graph here.
- `summary.json`, `REPORT.md`, four maps: structural counts, rain accounting,
  association coverage, censored observed branch spans, stage times and caveats.

Rain sums are equal-grid-cell sums over supplied RAINNC, not projection-area
corrected physical volume, and omit any unavailable RAINC contribution. Catalog
association fractions are not tracking accuracy; genesis, lysis, domain/time
boundaries and events legitimately affect them. Stage timings include different
amounts of work and follow a single warm-input schedule, so they are not a
controlled executor benchmark. C diagnostic time includes repeated seed checks
and a weak-component comparator, not just production envelope segmentation.

## Promotion gate

First inspect common-core retention, multicore grouping, weak-only exclusion,
centroid drift and event/lifetime changes on the paired 72 hours. Do not select
solely by retained rainfall, edge count or track longevity. Freeze the chosen
candidate definition before checking other dates; use independently chosen
windows, and later other years / future data with the same definition. The
current June window has already participated in tracking development and is
therefore a development sensitivity window, not independent validation.

0.1-mm/h envelopes are a separate next-stage experiment. This script deliberately
does not change their parameters or silently tune the frozen tracking algorithm.
