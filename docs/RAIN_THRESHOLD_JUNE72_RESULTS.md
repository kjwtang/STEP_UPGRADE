# 1.0 / 0.5-mm/h paired experiment: 1996 June 1, 72 hours

## Outcome: useful extent experiment, no threshold promotion

The fixed 1.0-mm/h tracking graph plus 0.5-mm/h extent (C) recovers most of the
additional rainfall represented by a direct 0.5-mm/h run (B), without changing
core identities or links. B substantially changes grouping, including large
objects. Keep A frozen; C is a promising **rainfall measurement layer**, not a
validated replacement for the storm definition. Neither B nor C is promoted.

The full 1749 x 2049 grid covers 72 hourly RAINNC intervals ending
1996-06-01 01:00 through 1996-06-04 00:00. Inputs were prepared previously from
73 cumulative snapshots with no negative increments. Timestamps follow the
explicit filename contract; internal time / latitude-longitude arrays were
not independently checked. This June window has participated in tracking
development: it is a sensitivity/development window, not independent truth.

## Rain accounting

| Quantity | A: direct 1.0 | B: direct 0.5 | C: 1.0 seeds + 0.5 partition |
|---|---:|---:|---:|
| Share of all supplied valid RAINNC rain | 79.750% | 88.465% | 87.733% |
| Labeled / assigned grid-time pixels | 6,646,063 | 10,661,017 | 10,279,916 |
| Extent relative to A | 1.000 | 1.604 | 1.547 |
| Object snapshots in the tracking graph | 16,769 | 17,935 | Exactly A |
| Graph edges | 5,526 | 6,050 | Exactly A |
| Split / merge / complex event records | 25 / 23 / 1 | 34 / 28 / 2 | Exactly A |
| Gap-continuation edges | 5 | 9 | Exactly A |
| Median observed branch span | 1 h | 1 h | Exactly A |
| Longest observed branch span | 56 h | 66 h | Exactly A |

C includes 99.173% of B's >=0.5-mm/h rain, while retaining 96.425% of its
eligible pixels. C excludes 381,101 weak eligible grid-time pixels and
240,355.96 mm-cells of weak rain: 0.732 percentage points of all provided rain.
That exclusion is small by rain contribution, but may still matter for storm
initiation/decay timing. C does not lengthen tracks by recovering that phase.

All sums use equal grid-cell weighting, not projection-area-corrected volume,
and supplied RAINNC only; unavailable RAINC is not included. The two thresholds
are rain-rate selection rules, not convective-core diagnoses.

## Grouping changes are substantial

B has 6,113 object snapshots with no A seed (34.08% of B objects). It also has
2,321 object snapshots containing more than one A seed. Those grouped snapshots
contain 7,268 A object snapshots, representing 84.69% of A's labeled rainfall;
730 of A's 1,006 >=1,000-cell snapshots participate. The largest B object
contains 24 distinct A seed labels. No A label was split across B labels in
this sample. These are grouping changes, **not adjudicated false merges**.

Only 126 seeded connected >=0.5-mm/h weak components contain multiple A seeds,
with at most four seeds per such component. In contrast, B's STEP radius-9
grouping joins 2,321 multi-seed objects. Thus the threshold change interacts
strongly with dilation-based joining; it is not just adding a weak-rain skirt.
The counts refer to different component definitions and do not provide a
one-to-one causal decomposition. Do not tune threshold and radius together
without a separate controlled experiment.

For 9,501 A object snapshots with bijective shared-core B correspondences,
2,359 projected links are retained, 753 are added, 138 are removed and five
retained links change event type. A total of 3,029 A edges and 2,938 B edges are
excluded from this comparison because one or both endpoints are ambiguous or
weak-only. This limited mapping is not whole-graph equivalence or accuracy.

Incoming labeled-rain association coverage is 88.09% in A and 89.96% in B;
one-observed-snapshot branches remain approximately 81.5% in both. Greater
coverage, longer tracks and more event edges are not proof of correctness.
Observed lifetimes are censored by the 72-hour window, domain boundaries and
intentional event-based branch termination.

## Expanded centroids must not silently feed the tracker

C's intensity-weighted core-to-envelope centroid drift is 0.424 cells at the
median, 2.585 at p90, 9.887 at p99 and 47.834 cells at the maximum. For cores
>=1,000 cells, p99 is 17.987 cells and the maximum is 35.924 cells. Drift is
therefore not exclusively a tiny-object effect. At nominal 4-km spacing, the
maximum all-object shift corresponds to approximately 191 km in grid space,
not a measured physical storm displacement.

The envelope/core cell-count ratio has median 1.75, p90 4.63, p99 28 and
maximum 1,680. The existing envelope is uncapped in distance; these extremes
need mask-level review. This is why C keeps the A graph and does not replace
its centroids, speed estimates, association scores or split/merge evidence.

## Reproduction and engineering checks

- A reproduces the prior frozen June-1 reference exactly: every identified and
  tracked frame hash, complete object catalog, full edge scores/overlaps,
  events and final family roots. A separate final audit also matches saved
  timestamps, physical grid attributes, hourly cadence and units.
- Both experimental graphs pass node/branch uniqueness, non-reuse, forward
  time/gap and event/checkpoint-reference checks.
- Every frame passes exact A-seed equality, partition core ownership,
  partition/system assigned-support equality, dry/missing barrier checks and
  `C assigned rain + unassigned eligible rain = B rain` conservation.
- Synthetic tests cover weak bridges, weak-only rain, unequal numeric labels,
  many-to-one exclusions, changing chunk boundaries, dry/NaN windows,
  negative-rate rejection and reference timestamp mismatches.

The complete paired diagnostic took 546.15 seconds locally. A identification /
tracking / output stages took 38.74 / 53.26 / 13.83 seconds; B took
35.10 / 78.62 / 21.07 seconds. C extent and accounting diagnostics took
256.38 seconds. C timing includes repeated seed reproduction, weak-component
comparison and geometry/accounting; it is **not** a production envelope cost.
The A/B runs used eight frame workers in six-frame batches, hence no more than
six useful frame tasks per batch. Testing overlapped part of this run, and
plotting initialized a font cache; timings do not establish controlled speed
ratios or RCC scaling.

Local results are under `results/threshold_paired_june72_run1/`: four preview
maps, per-frame accounting, shared-pixel membership, weak-component membership,
envelope geometry, projected edge changes and complete A/B stream catalogs.
`configuration.json` retains the actual run's script/input hashes; the first
real run preceded the addition of the reusable time/grid audit helper. The
helper was then applied to the saved real catalogs and passed all eight checks.
See [the runnable protocol](RAIN_THRESHOLD_SENSITIVITY.md) for reproduction.
No raw inputs, generated arrays, catalogs or images are committed to GitHub.

## Next step, not performed here

Keep all tracking parameters fixed. Check the same A/B/C definitions on dates
not used to choose the threshold, including organized bands, isolated cells
and weak-rain periods. Inspect the largest grouping and extent/drift cases;
do not select a winner solely by rain contribution. A minimum stronger seed,
seed area/duration rule, distance cap, radius change and 0.1-mm/h extension
must each be explicitly named separate sensitivity experiments. Rain-only
tests cannot certify convection; independent years/future members and stronger
physical evidence remain necessary for scientific promotion.
