# Chang et al. (2016) S1 reference implementation

`step/chang2016.py` implements the identification and tracking algorithm of
Supplemental Material S1 of Chang, Stein, Wang, Kotamarthi and Moyer (2016),
J. Climate 29, 8355-8376 (doi:10.1175/JCLI-D-15-0844.1, supplement
10.1175/JCLI-D-15-0844.s1). It is additive: `step.identification`,
`step.tracking`, the reference preset and the streaming runner are unchanged.
Use it as the paper baseline against which STEP_UPGRADE choices are measured.

## Source facts (supplement S1, manuscript sections 2-4)

| Item | Paper value |
|---|---|
| Sampling | 3-hourly; Stage IV hourly data aggregated to 3 h |
| Wet cutoff | > 0.1 mm per 3 h (about 0.033 mm/h) |
| Structuring sets | disks, radius 60 km (JJA) and 120 km (DJF) |
| Large region | at least one cell left after erosion |
| Large grouping | opening, then almost-connected labelling of the opened regions |
| Small regions | each dilated individually; joins the large segment sharing most cells; the rest are grouped by almost-connected labelling |
| Similarity | sum_s sum_u w_i(s) w_j(u) exp(-phi d(s,u)), w = rain / segment rain sum |
| phi | 1/120 km^-1 (d is Euclidean distance in the model/Stage IV Cartesian grid) |
| tau | 0.05 (JJA), 0.01 (DJF); link requires sim > tau |
| Motion rule | centroid displacement < 120 km, or angle to the previous displacement < 120 degrees |
| Choice | qualifying previous segment with the largest number of grid cells |
| Events | one predecessor per segment; several segments may share a predecessor (split) |

Units: intensity is mm/h (mean over the segment's cells; size-weighted over
events for the aggregate factor); size is cells x 144 km2 (12 km) or
x 16 km2 (4 km). The paper never uses mm/h/km2. The product
Intensity x Size x Duration x Number is a rain volume (mm km2;
1 mm km2 = 1000 m3). `factorize()` returns the four factors and checks this
identity.

## Parameters per grid

| Quantity | 12 km grid | 4 km grid |
|---|---:|---:|
| 60 km disk radius | 5 cells | 15 cells |
| 120 km disk radius | 10 cells | 30 cells |
| phi per cell | 0.1 | 1/30 |
| 120 km displacement | 10 cells | 30 cells |

For comparison, the original Python STEP tutorial uses `phi=0.003` per cell
and `tau=0.7`; the STEP_UPGRADE reference preset uses a 9-cell (36 km at
4 km) bridging radius and a 1 mm/h cutoff. None of these equal the paper.

The paper's tau, phi and 120 km rule were tuned on 3-hourly fields. With
hourly fields a storm moves about a third as far per step, so similarities
are higher and the 120 km rule is three times looser per hour. The presets
`configs/chang2016_s1_{summer,winter}_3h.json` therefore average hourly
input to 3-hourly rates. Running the paper parameters on hourly frames is an
extrapolation and should be reported as such.

## Choices the supplement leaves open

* `strict=True`: wet means rain > threshold, as worded in the paper.
* `erosion_border='ignore'`: cells outside the domain do not erode regions,
  matching `skimage.morphology.erosion` in the original Python STEP.
* A large region whose opened cells fall into several almost-connected
  clusters joins those clusters, because the region itself is contiguous.
* Small-region assignment ties go to the lowest segment id; choice ties go to
  higher similarity, then the earlier previous segment.
* `no_history='reject'`: when the displacement is >= 120 km and the previous
  segment has no predecessor, the angle cannot be computed and the link is
  refused. The original Python STEP accepts such links at the second time
  step (`vec_angle = 1`); use `no_history='accept'` to reproduce that.
  A zero previous displacement also makes the angle undefined (refused).
* `selection='largest'` follows S1.2. The manuscript text speaks of the
  "most morphologically similar" event; `selection='most_similar'` is a
  sensitivity option.
* Missing data (NaN or `valid_mask=False`) is dry and blocks almost-connected
  dilation; the paper had no missing data.
* Supplement S2 writes x = cos(lat) sin(lon) and y = cos(lat) sin(lon); the
  first is a typo for cos(lat) cos(lon). Tracking uses grid centroids as S1.2
  specifies, so S2 does not enter here.

## Differences from the original Python STEP (RDCEP/STEP ef75c08)

* Original identification erodes and dilates integer label images, so results
  depend on label values; this module uses binary masks.
* Original small-region assignment dilates groups of small regions; S1.1
  dilates each small region individually.
* Original tracking compares a current segment with the union of all
  previous-step cells sharing an event id; S1.2 compares segments.
* Original `vec_angle` can carry over from an earlier candidate when the
  previous storm is absent at t-2.
* Original similarity allocates several dense U x U matrices (U = union of
  cells); see docs/ORIGINAL_STEP_ANCHOR.md for the RCC memory failure.

## Differences from the STEP_UPGRADE lineage tracker

The lineage tracker scores IoU, centroid distance and an intensity ratio,
solves one-to-one assignments, and creates new branches at splits and merges.
Neither `branch_id` nor `family_id` equals a paper event: branches end at
splits (too short), families union merged events (too long and too few).
Size, duration and number factors depend on that definition, so compare
factors only with a stated event definition.

## Computation

* Disk erosion and dilation use exact Euclidean distance transforms; cost does
  not grow with the radius. Tests check equality with scipy binary morphology.
* Similarity is exact. Pairs whose bounding-box gap g satisfies
  exp(-phi g) <= tau cannot exceed tau and are skipped (exact bound; for
  tau = 0.05 and phi = 1/120 km^-1 this is g >= 359.5 km). Small pairs are
  summed directly in blocks; larger pairs use a linear FFT convolution with the
  full, untruncated kernel on the field of the larger segment, shared by all
  its partners. Kernel spectra are cached. Memory is linear in domain size.
* Synthetic 1010 x 1634 fields (single core, this test environment): about
  0.6 s per frame identification; tracking about 4.4 s per frame with
  ~100 segments (60 km disk) and 0.4 s per frame with a few large segments
  (120 km disk). These are synthetic timings, not production estimates.
  `fft_workers` enables threaded FFTs.

## Use

```bash
python scripts/run_chang2016_reference.py rain_mm_h.npy \
  --preset configs/chang2016_s1_summer_3h.json --grid-km 4 \
  --output-dir results/chang2016_summer_01
```

Input is a (time, y, x) array of hourly rain rates in mm/h on a uniform grid.
Outputs: `segment_labels.npy`, `event_labels.npy`, `segments.csv`, and
`summary.json` with the four factors and their units. `--cell-area-km2`
overrides the default `grid_km**2` (the paper's constant cell area).

`step/identification_fast.py` is a separate, opt-in accelerator for the
existing STEP_UPGRADE identification with disk footprints; it returns
identical labels (tests/test_identification_fast.py).
