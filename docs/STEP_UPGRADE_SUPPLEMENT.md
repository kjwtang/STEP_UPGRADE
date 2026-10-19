# STEP_UPGRADE algorithm notes (draft supplement)

Status: working draft, 2026-10-18. Informal; written in the spirit of the
Chang et al. (2016) supplement S1 so that the hourly 4 km algorithm can be
described, reviewed and cited later. Numbers come from local experiments on
CSTM WRF output (4 km, hourly, 1749 x 2049 cells) and are illustrative, not
final results.

Everything described here is opt-in. The `v0.4.0b3` modules
(`step/identification.py`, `step/tracking.py`), the reference preset
`configs/rainfall_lineage_4km_1h_reference_v1.json` and their defaults are
unchanged. With default arguments, every new class or function reproduces the
existing behaviour exactly; tests check this.

Related: `docs/CHANG2016_REFERENCE.md` (the paper algorithm and its units).

## 1. Why these changes

An audit of the preset on real data found these problems, in order of
impact:

1. **Chaining in identification.** The preset dilates the whole 1 mm/h mask
   by a 9-cell (36 km) disk and labels the result. Any chain of cells closer
   than 36 km becomes one object. In 72 h of June 1996, 575 objects had a
   bounding-box diagonal > 1000 km, and they held 65.8% of the rain. The
   largest object spanned 4828 km and consisted of 1283 disconnected pieces.
2. **IoU scoring with one-to-one assignment.** The overlap term of the score
   is IoU, which is tiny for pairs of very different size. For example, a
   50-cell child lying entirely inside a 3000-cell parent has IoU of about
   0.02. Once the parent is assigned to its main successor, such children can
   only start new branches. The split rule needs an intersection of at least
   15% of the parent, so small fragments are not splits either.
3. **Flicker.** Small, fast-changing cells that overlap only slightly
   (median 4 cells, IoU 0.04) score below tau = 0.35. The gap rule
   (`gap_tau` 0.45, `gap_min_pixels` 16) almost never fires: 11 gap links in
   72 h.
4. **Event definition.** `branch_id` ends at every split or merge, and
   `family_id` unites everything that ever split or merged. On the test data
   this turned out to matter little (see section 4).
5. **Grid- and step-dependent parameters.** Radii, displacement and pixel
   minima are in cells and frames. Some 16-cell minima are hard-coded in
   `tracking.py` (lines 552, 588, 625 and 754) and cannot be changed through
   the configuration.
6. **The 1 mm/h threshold removes 16-25% of the rain**, and the removed
   share changes between climates (for example 19.5% -> 16.2% in August,
   1996 -> 2046).

Items 1-4 have opt-in fixes, described below. Items 5 and 6 are open.

## 2. Identification: `step/identification_regions.py`

`identify_regions(data, radius_cells, threshold=1.0, ...)` follows the
four-step structure of S1, with two generalizations.

1. Contiguous wet regions (8-connectivity), `rain >= threshold`.
2. A region is **large** if it keeps at least one cell after erosion by a
   disk of `radius_cells` (the *core radius*).
3. Large regions are grouped (`large_join='almost'`, as in S1): their
   openings are dilated by the *bridge radius* and connected.
   `large_join='region'` keeps every large region separate.
4. Each small region is dilated individually by the bridge radius and joins
   the large segment it covers most (ties: lowest id). A small region never
   connects two large segments.
5. Small regions that reach no large segment are almost-connected among
   themselves (`small_join='almost'`, as in S1) or kept separate
   (`small_join='region'`).

Options:

| Argument | Default | Meaning |
|---|---|---|
| `radius_cells` | required | core radius used by the erosion test (step 2) |
| `bridge_radius_cells` | `radius_cells` | dilation radius for steps 3-5; S1 uses one radius for both |
| `large_join`, `small_join` | `'almost'` | S1 joining steps, or `'region'` |
| `edge_cells` | 0 | band along the domain border treated as invalid; no rain is labelled there and no bridge crosses it |
| `min_size` | 1 | drop objects smaller than this |

`edge_mask(shape, n)` returns the same band as a boolean mask. It can be
passed as `valid_mask` to `identify_disk` to apply the band to the preset.

With defaults, the output equals `chang2016.identify_frame(strict=False)`
pixel for pixel.

### What the experiments showed (72 h, 1 mm/h)

| Identification | Objects | Objects > 1000 km | Rain in them | Longest (km) |
|---|---|---|---|---|
| preset, disk 9 | 14,268 | 575 | 65.8% | 4828 |
| preset, disk 5 | 29,002 | 358 | 42.0% | 4029 |
| regions, core 9 / bridge 9 (plain S1) | 16,559 | 484 | 40.3% | 4316 |
| regions, core 9 / bridge 9, `small_join='region'` | 135,882 | 148 | 16.6% | 2009 |
| **regions, core 3 / bridge 9** | **24,088** | **287** | **25.0%** | **2685** |
| regions, core 3 / bridge 5 | 41,043 | 168 | 18.3% | 2095 |

* With plain S1, most of the remaining chaining comes from step 5, not from
  step 3: switching off `large_join` changes almost nothing (484 -> 483).
  At 1 mm/h few regions are wide enough to contain a 36 km disk, so most
  rain is classified as small and is then bridged by step 5.
* Switching off step 5 removes the chaining but fragments the field (about
  10 times more objects).
* A smaller core radius with the same bridge radius is the compromise.
  Changing the method is much more effective than shrinking the preset disk.
  For a similar object count, the regions method leaves far less rain in
  chained objects.
* Visual checks (case maps) confirm that the preset's longest objects are
  scattered convective cells joined into one chain. With core 3 / bridge 9,
  the 4828 km object at 1996-06-17 19 UTC splits into 58 objects.
* Cost: about 0.45 s per frame, against 0.15 s per frame for the preset
  (single process, full domain). A KD-tree replacement for the per-region
  attachment loop gave identical labels but was slower, so it was dropped.
  The time is spent in full-domain distance transforms.

Domain edge: mean rain within 5-30 cells of the border is 1.6-2.2 times the
interior mean. Straight, one-cell-wide strips along the east and south edges
appear in the case maps. This looks like the WRF lateral relaxation zone
(partly ocean). The edge band holds 1-5% of the rain; masking 15 cells
changes preset N by less than 1%.

## 3. Tracking: `step/tracking_coverage.py`

`CoverageTracker` subclasses `Tracker`, and `track_with_graph_coverage`
mirrors `track_with_graph`. The new options are independent and can be
combined.

| Option | Effect |
|---|---|
| `overlap_metric='iou'` (default) | unchanged |
| `overlap_metric='child_coverage'` | the overlap term of the adjacent score uses intersection / child area instead of IoU |
| `overlap_metric='coverage'` | uses intersection / min(area) instead; not recommended, see below |
| `share_parent=True` | after the normal update, a child with no edge and at least `share_min_coverage` (0.5) of its cells inside one previous object joins that object's branch (edge `share`); a branch can then hold several nodes at one time, as S1 allows |
| `rescue_overlap=True` | pairs unlinked children with parents that have no successor when each is the other's largest raw overlap (at least `rescue_min_pixels` cells): first adjacent (`overlap_continue`), then one frame apart (`gap_overlap`); stays one-to-one |

The split/merge rules and the strong- and ranked-overlap rules are not
changed. Children joined by `share` or `rescue` carry no velocity, because
their centroid is not a motion estimate of the parent. Checkpoints record
the options, and resuming with different options raises an error.

### What the experiments showed (72 h)

"Born inside" counts nodes with no edge and at least 50% of their cells
inside a previous object. "Wrong parent" counts continue edges whose parent
is not the largest overlap of the child.

| Segments | Setting | Branches | Duration (h) | Born inside | Wrong parent |
|---|---|---|---|---|---|
| preset | default | 9232 | 1.55 | 353 | 4 |
| preset | share | 8884 | 1.56 | 5 | 4 |
| preset | coverage | 8056 | 1.77 | 207 | 141 |
| preset | child_coverage | 8490 | 1.68 | 211 | 6 |
| preset | child_coverage + share | 8283 | 1.68 | 4 | 6 |
| paper (S1, 60 km) | default | 5391 | 1.79 | 839 | 21 |
| paper (S1, 60 km) | child_coverage + share | 4305 | 2.03 | 86 | 32 |

* `share_parent` fixes item 2 and leaves every existing edge untouched.
* `coverage` hands large children to small parents: when a small parent
  lies inside a large child, the overlap coefficient is also 1. Wrong parents
  rise from 4 to 141, so this option should not be used. `child_coverage`
  does not have this problem.
* `rescue_overlap` links many more cells (N -6 to -14%), but its evidence is
  thin. The median overlap is 3 cells (7% of the child), and the 90th
  percentile of centroid displacement is 200-400 km/h. Treat it as a
  sensitivity option only.

## 4. Events: `step/lineage_events.py`

`paper_event_ids(graph)` assigns S1 events from a lineage graph:

* a node without incoming edges starts an event;
* a single incoming edge passes the parent's event on, so split children
  share it;
* a merge keeps the event of the parent with the most cells (ties: lowest
  node id).

`event_factors(...)` returns N, the sum of event lifetimes, size and
duration.

Finding: on every tested setting, the number of paper events equals the
number of families. All merges re-join objects that had split from a common
ancestor, so each family has exactly one root. The roughly 2x gap between
lineage tracking and the paper rules on identical segments (4094 against
2099 events in 72 h) therefore comes from *fewer links*, not from the event
definition. Part of that gap is the paper rule's own over-linking: 12-16%
of its hourly links move the centroid by 300 km or more.

## 5. Hourly variant of the paper tracker: `step/chang2016_hourly.py`

`track_hourly` keeps S1 events and similarity, and replaces three rules:
mask-to-mask distance (default 40 km) instead of the centroid and angle rule,
the most similar parent instead of the largest, and an optional one-frame
gap. Compared with the paper rule using most-similar selection, its N and
duration stayed within a few percent in the summer windows. In the two
winter windows it gave 20-24% more events per day. Its main benefit is that
every link is spatially adjacent by construction. It is 2-3 times slower
than the paper tracker, because of per-pair KD trees.

## 6. Recommended configuration

```python
from step.identification_regions import identify_regions
from step.tracking_coverage import track_with_graph_coverage
from step.lineage_events import paper_event_ids

labels = identify_regions(rain_mm_h, 3, threshold=1.0,
                          bridge_radius_cells=9, edge_cells=15)
tracking = dict(preset["tracking"])            # the reference preset
tracking["km"] = tracking.pop("max_displacement")
raster, graph = track_with_graph_coverage(labels, rain_mm_h,
                                          overlap_metric="child_coverage",
                                          share_parent=True, **tracking)
events = paper_event_ids(graph)                # optional, S1-style ids
```

Results over 8 windows: summer weeks of June, July and August in 1996 and
2046, and two 3-day DJF windows in 1996.

| | Preset | Recommended |
|---|---|---|
| Summer rain in objects > 1000 km | 65-82% | 20-40% |
| Winter, same | 89-91% | 57-65% (mostly real stratiform systems) |
| Summer branches per day | 2940-4305 | 4338-5842 |
| Summer size (10^4 km2) | 0.56-1.00 | 0.37-0.59 |
| Summer duration (h) | 1.39-1.56 | 1.59-1.74 |

Chaining biases the change signal. In July, the size change 1996 -> 2046 is
+37.6% with the preset and +26.1% with the recommended setting; the change in
N is -14.4% and -5.2%. June and August differ less. Each window is one week
of one year, so these are method comparisons, not climate results.

## 7. Open items

* Item 5: parameters in km and hours, and the hard-coded 16-cell minima in
  `tracking.py`.
* Item 6: the threshold. For the main analysis, a low threshold (about
  0.033 mm/h) with a 1 mm/h intense-rain subset is suggested.
* Identification speed: restrict distance transforms to wet bounding boxes.
* Choose the core radius on physical grounds. Core 3 / bridge 9 keeps larger
  organized systems; bridge 5 gives cleaner cell clusters.
* Known failing test: `test_split_shares_event_and_merge_continues_largest`
  in `tests/test_chang2016.py`. The case sits exactly on the 120 km threshold
  (120.00000000000003 km, which is not < 120), so it is a floating-point
  boundary in the test, not a tracker error.
* Data note: in the 2046-07-21 files used here, eight hours came from a
  different model run (by `START_DATE`), which breaks cumulative RAINNC.
  Check for run overlaps before differencing.
