# 0.5 / 0.1-mm/h extents: frozen June 72-hour tracking

## Outcome

0.1-mm/h extent is useful for rainfall coverage but is **not promoted** as a
storm definition. The unchanged 1-mm/h tracking graph is reused for both extents.
Large, weak extents attached to tiny seeds require explicit scientific treatment;
rainfall recovery alone does not resolve their identity. No tracking score, seed
threshold, radius, event rule, stronger-core criterion or distance cap was tuned.

This is the same full 1749 x 2049 grid and hourly RAINNC window ending
1996-06-01 01:00 through 1996-06-04 00:00 used in earlier development. It is a
paired sensitivity test, not an independent date/year validation.

## Rain and extent accounting

| Quantity | 0.5 extent | 0.1 extent |
|---|---:|---:|
| Assigned fraction of all provided rain | 87.733% | 95.552% |
| All threshold-eligible rain fraction, including unseeded components | 88.465% | 97.149% |
| Assigned grid-time pixels | 10,279,916 | 20,106,176 |
| Unassigned eligible grid-time pixels | 381,101 | 2,375,417 |
| Multi-seed wet-component snapshots | 126 | 558 |
| Maximum distinct seed labels in one wet component | 4 | 16 |
| Extent/core cell-count ratio: median | 1.75 | 3.91 |
| Extent/core cell-count ratio: p99 | 28.0 | 325.67 |
| Extent/core cell-count ratio: maximum | 1,680 | 11,848 |
| Rain-weighted core-to-extent centroid drift: median, cells | 0.424 | 0.758 |
| Centroid drift: p99, cells | 9.887 | 28.902 |
| Centroid drift: maximum, cells | 47.834 | 132.273 |

The lower extent adds 7.819 percentage points of supplied rainfall and multiplies
assigned pixels by 1.956. All rain sums use equal grid-cell weights and provided
RAINNC, not projection-corrected physical volume or an unavailable RAINC field.

The extra support contains 9,579,541 pixels below 0.5 mm/h **and 246,719 existing
>=0.5-mm/h pixels** that were previously coreless/unassigned but become reachable
through weaker bridges. Those existing stronger pixels contribute 5.95% of the
added rain. Thus the change is not only a skirt of lower-intensity pixels.

Previously assigned pixel ownership changes: **zero** in all 72 frames. Exact
seed identity and nested assigned support both pass. This measured result should
not be generalized to different segmentation fields, seed choices or algorithms.

## Extreme cases need mask review

| Frame / interval end | Seed / branch | Seed cells | 0.1 extent cells | Drift, cells | Why inspect |
|---|---|---:|---:|---:|---|
| 28 / June 2 05:00 | local 229 / branch 4902 | 1 | 11,848 | 72.685 | Maximum area ratio; tiny seed determines a broad weak extent |
| 10 / June 1 11:00 | local 260 / branch 217 | 193 | 15,847 | 132.273 | Maximum rain-weighted centroid drift |
| 4 / June 1 05:00 | local 244 / branch 573 | 41 | 18,439 | 93.352 | Farthest assigned pixels from any seed: 275.363 cells |

At nominal 4-km spacing, the maximum centroid drift corresponds to about 529 km
in grid space. It is **not** storm motion. The 275-cell diagnostic is distance to
the nearest pixel of ANY retained seed, not its assigned seed and not a wet-path
length. Do not present either as a measured geographic distance or use it to set
a meteorological distance cap automatically. Full-grid maps preserve grid axes;
latitude/longitude placement has not been inferred.

Extreme growth is not exclusively a single-pixel issue. Among the 1,006 seed
snapshots with >=1,000 cells, 0.1 extent/core ratio has median 2.35, p99 10.87 and
maximum 19.57; centroid drift has p99 32.50 and maximum 61.97 cells. The 1,000-cell
bin is a descriptive analysis subset, not a new minimum storm area.

These cases do not establish false links or false extent assignments. They show
why "1-mm/h seed + unrestricted 0.1-mm/h extent" needs an explicit measurement
interpretation and why expanded centroids must not drive frozen tracking.

## Reproducibility and tests

- All 72 reconstructed seed raster hashes match saved A exactly; chronology,
  grid, units, catalog seed IDs and areas agree.
- Saved graph/catalog/checkpoint sources remain unchanged. Both extents reference
  the same 16,769 nodes, 5,526 edges and 11,342 branch IDs. The tracker was not rerun.
- 0.5 accounting and object geometry reproduce the earlier C experiment. Its
  original C raster hashes were not saved, so unavailable pixelwise equivalence
  is not claimed.
- Seed preservation, wet/NaN barriers, seeded-component support, nested support
  and rain conservation pass in every frame.
- The implementation full suite passed 152 tests with one skipped. Existing
  dependency warnings (NumPy/netCDF4/xarray) were emitted; no test failed.
- New tests cover noncontiguous seed values, immutable seeds, connectivity and
  optional cap compatibility, chunk invariance, coreless rain, invalid input,
  checksum/hash mismatch and failed-run completion guards.
- An additional bounded 2,000-case synthetic ownership search found no existing
  pixel reassignment; this is a diagnostic observation, not a proof for all cases.

The actual local run took 275.26 s including startup, checks and plots. Seed
replay took 43.22 s; both extents plus accounting/geometry took 165.45 s. Tests
overlapped part of the run and plotting initialized a temporary font cache, so
these timings are not a controlled speed comparison. Seed replay was parallel;
watershed/statistics were serial per frame. No 16/32-worker scaling claim is made.

Results: `results/frozen_extents_june72_run1/`, including seven full-grid maps,
immutable frame/node/branch mappings, ownership/support accounting and ranked
extremes. Raw data, arrays, catalogs and generated figures remain Git-ignored.
The saved configuration fingerprints the exact experiment code. Afterwards only
ownership-warning wording and reference-manifest coverage were tightened; the
segmentation/statistical calculations were not changed. See
[the protocol](FROZEN_EXTENT_SENSITIVITY.md) for reproducible commands and scope.

## Decision and next experiment

Keep tracking frozen and retain 0.5 as the existing experimental measurement
comparison, **not** a newly certified default. Keep 0.1 as a separately named
diagnostic extent. Do not add an arbitrary distance cap or discard small-seed
cases merely because they look extreme.

Next choose an independent date with broad weak rain and another with isolated
strong objects. Assess the same fixed definitions, inspect the extreme masks and
report unseeded onset/decay rain. If a stronger-core interpretation is required,
compare seed area, rate and lifecycle strength explicitly and account for excluded
rain. An optional coherent weak-system grouping and a per-core ownership layer
should remain separate; neither is automatically the original STEP event.
