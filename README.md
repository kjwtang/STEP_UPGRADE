# STEP_UPGRADE

STEP_UPGRADE turns gridded rainfall into a time-continuous catalog of rain objects:
identify each hourly object, track its movement and relationships, and retain
continuations, splits, merges, and short gaps across files, chunks, and months.

Group beta: **v0.4.0b3**. The recommended workflow uses hourly 4-km data and the
[1-mm/h reference preset](https://github.com/kjwtang/STEP_UPGRADE/blob/v0.4.0b3/configs/rainfall_lineage_4km_1h_reference_v1.json).
This README describes that beta. The program operates on rainfall, not cloud imagery.

## 1. Input: cumulative RAINNC directly

**The recommended runner accepts cumulative `RAINNC` directly. You do not need
to prepare single-hour rainfall or merge files into one large NetCDF first.**
It accepts an hourly snapshot directory, a single multi-time NetCDF, multiple
multi-time files, or a directory containing those files. Multi-time fields may
be cumulative amounts, hourly interval amounts, or rain rates; the reader
detects their type from the variable name, units, and metadata.

It reads consecutive snapshots and computes hourly rainfall internally:

```text
rain_rate[t] = (RAINNC[t] - RAINNC[t - 1]) / 1 hour
```

For 24 rainfall intervals, supply 25 consecutive snapshots. The first snapshot
is the predecessor needed for differencing; it is not an extra tracked hour.
Each interval is timestamped by its ending snapshot.

| Input item | Required format |
|---|---|
| Hourly snapshot files | Actual NetCDF named `cstm_d01_YYYY-MM-DD_HH:MM:SS.nc` or `.rain`; cumulative `RAINNC` in mm |
| Snapshot dimensions | `(south_north, west_east)`, optionally preceded by singleton `Time` |
| Multi-time files | NetCDF precipitation with three axes, normally `(time,y,x)` or `(Time,south_north,west_east)`; axis order may differ |
| Multi-time variable | `RAINNC` selected automatically; otherwise one identifiable precipitation variable, or `--variable NAME` |
| Time | Default interval is 1 hour; decoded calendar times must be consecutive; index-only time uses the hourly assumption below |
| Grid | Same shape/projection in every file; the reference preset requires `DX=DY=4000` meters |
| Coordinates | `XLAT/XLONG` arrays are not required by the reference runner |
| Sequence | One climate/member sequence per run; do not mix members or duplicate timestamps |

Hourly CSTM snapshot directories require global grid attributes:
`DX, DY, MAP_PROJ, CEN_LAT, CEN_LON, TRUELAT1, TRUELAT2, STAND_LON,
MOAD_CEN_LAT, POLE_LAT, POLE_LON`. Use actual source metadata. Multi-time files
require `DX`/`DY`; remaining projection attributes are retained and compared when
present. If spacing is absent, supply the actual spacing with `--grid-km 4`.
This explicitly declares spacing; it does not generate latitude/longitude or a projection.

The runner scans the input directory recursively. Missing hours, duplicate
timestamps, changing grids, negative accumulation, or cumulative decreases/resets
stop the run. Resolve those inputs before rerunning; they are not silently
zero-filled. Only the selected variable is processed; RAINC is not added automatically.

### Automatic type detection and time assumptions

- `RAINNC` with amount units is treated as cumulative. Explicit `rain_kind`
  attributes and cumulative descriptions are also recognized.
- Rate units such as `mm/h` or `kg m-2 s-1` identify a rate and are converted
  to mm/h. Amounts marked `rain_kind=interval` or `cell_methods="time: sum"`
  are hourly interval amounts. No differencing is applied to interval/rate input.
- An ambiguous amount field is not guessed from its numerical values. Use
  `--rain-kind cumulative`, `interval`, or `rate`. Already differenced data still
  named `RAINNC` must declare `rain_kind=interval` or use `--rain-kind interval`.
- Numeric `time=0,1,2,...`, or a missing time coordinate, defaults to hourly
  indices. Numeric spacing must be hourly; non-hourly calendar data is rejected.
  Units are required; use `--units mm` only if this is the actual source unit.
- Calendar-timed files are ordered by their internal times. Index-only files
  use explicit argument order, or lexicographic path order in a directory.
  Zero-based indices may restart in each file. Such file boundaries are assumed
  continuous and recorded as assumptions, not verified dates. Use zero-padded,
  chronologically ordered filenames or pass the files explicitly in order.
- For index-only input, `--time-origin 1996-01-01T00:00:00` gives the first global
  sample a calendar date. Without it, `timestamps` is null and `sample_indices`
  identifies samples; the program does not invent calendar dates.

For cumulative input, differencing crosses file boundaries. Only the first
sample of the entire planned sequence is a predecessor rather than a tracked
hour. For example, 741 samples produce 740 hourly rainfall frames. If accumulation
resets at a monthly/restart boundary, correct that reset upstream; the first
sample in the next file is not silently substituted for its hourly rainfall.
Duplicate calendar snapshots are rejected; supply one copy of each time.

**Already have hourly rainfall?** Multi-time NetCDF interval/rate input can use
this same runner. The Python functions in section 5 take prepared mm/h arrays.

## 2. Install and check the environment

Use a separate checkout and virtual environment. Original STEP and STEP_UPGRADE
both import as `step`, so do not install them in the same environment.
Recommended Python: **3.13.5**; minimum: **3.12**. The runner uses Linux/macOS
POSIX process and file-lock facilities.

```bash
git clone --branch v0.4.0b3 --single-branch \
  https://github.com/kjwtang/STEP_UPGRADE.git STEP_UPGRADE_beta_040b3
cd STEP_UPGRADE_beta_040b3
python3.13 -m venv .venv-step-beta
source .venv-step-beta/bin/activate
python -m pip install --only-binary=:all: -r requirements-beta.txt
python -m pip install --no-build-isolation -e .
python -m pip check
python -c "import step; from importlib.metadata import version; print(version('step-upgrade')); print(step.__file__)"
```

Expect version `0.4.0b3` and an import path in this checkout. Keep the checkout:
the installation above is editable. Do not change source or dependencies while
a run is active or waiting to resume. Pins cover top-level dependencies; save
`pip freeze` with your run.

Check installation without research data:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -u scripts/beta_smoke.py --output-dir results/beta_smoke_01 --workers 2
```

Expected output: 12 frames, 11 nodes, 10 edges, one gap edge, ten equality
checks set to true, and `results/beta_smoke_01/SUCCESS`.
For RCC module/conda setup, see the
[installation guide](https://github.com/kjwtang/STEP_UPGRADE/blob/v0.4.0b3/docs/GROUP_BETA_0_4_0b3.md).

## 3. Main workflow: identify and track cumulative rainfall

**Entry point:** `scripts/run_reference_stream.py`.

**Use:** pass an hourly snapshot directory or one/multiple multi-time files. The command reads
and differences rainfall, identifies objects, advances tracking in time order,
and publishes tables and checkpoints per chunk.

```bash
python -u scripts/run_reference_stream.py /PATH/TO/ONE_SEQUENCE \
  --preset configs/rainfall_lineage_4km_1h_reference_v1.json \
  --start 0 --hours 24 --sequence-id member0_1996 \
  --workers 4 --chunk-frames 6 --id-executor fresh \
  --save-labels \
  --output-dir results/rain24_01
```

Inspect a monthly file before processing (no output directory is created):

```bash
python scripts/run_reference_stream.py /PATH/TO/month.nc --inspect
```

Process all available hours across multiple files in one continuous execution:

```bash
python -u scripts/run_reference_stream.py /PATH/TO/june.nc /PATH/TO/july.nc \
  --sequence-id member0_1996 --workers 4 --chunk-frames 6 \
  --output-dir results/multiple_months_01
```

The same command accepts just `/PATH/TO/month.nc` or a directory of monthly files.
Omit `--hours` to process the remaining sequence; set it for a bounded trial.
Add `--grid-km 4` only if DX/DY are absent and the data really has 4-km spacing.

| Option | Meaning / usage |
|---|---|
| `--start` | Global sample index: cumulative predecessor, or first interval/rate frame |
| `--hours` | Number of hourly frames; omitted means all remaining available intervals |
| `--inspect` | Print detected variable, kind, units, axes, timing assumptions and grid; no tracking |
| `--variable` / `--rain-kind` | Resolve variable selection or ambiguous amount metadata explicitly |
| `--dims TIME Y X` | Declare three axis names when automatic axis detection is insufficient |
| `--units` | Explicit source units; conversion to mm/h is recorded in chunk metadata |
| `--time-coordinate` / `--time-origin` | Select a 1D time variable, or date the first global index-only sample |
| `--grid-km` | Declare missing multi-time grid spacing; cannot contradict source DX/DY |
| `--sequence-id` | Namespace for one continuous climate/member sequence and its IDs |
| `--workers` | Parallel identification workers; tracking state advances sequentially |
| `--chunk-frames` | Rainfall frames processed per chunk; begin with 6-24 |
| `--save-labels` | Save per-chunk identification/branch arrays; omit for less disk I/O |
| `--output-dir` | New directory for the run; existing results are not overwritten |
| `--no-progress` | Disable the hourly percentage/ETA display |
| `--stop-after-chunks` | Pause after a specified number of committed chunks |

The preset uses a >=1-mm/h wet mask, a 9-cell morphological bridging radius,
and minimum object size of one wet cell. Grouping assigns nearby wet fragments
to a label but does not add artificial rain pixels to its output. A label may
contain disconnected fragments. The 9-cell radius is not a maximum object
separation. Missing pixels act as barriers.

Tracking retains continuation, split, merge, complex many-to-many events, and
gap continuation. Gap recovery can bridge one frame where an object is absent;
both endpoints must have at least 16 wet cells under the reference preset.
It cannot replace a missing raw cumulative file.

**Output:** CSV object/relationship tables, JSON metadata/checkpoints, optional
NumPy label arrays, and a final summary. See the exact schemas in section 4.
The runner does not generate preview plots automatically.

### Resume and cross-month continuity

Repeat the same command and append `--resume` after an interruption. Keep
the input plan, `--hours`, preset, sequence ID, source, environment, and
`--save-labels` setting unchanged. Worker count and chunk size may change.

Include adjacent months in one planned input sequence; carry tracking state
rather than tracking each month independently and concatenating its IDs.
A new climate/member needs a new sequence/output directory.
The input file boundaries and processing chunk boundaries need not coincide.
For cumulative input, each chunk rereads its predecessor, including across files;
tracking continues from the committed checkpoint without renumbering.
Resume requires the same planned files and interval range. This is not a watch
mode: newly appended hours require a separately planned workflow, not changing
`--hours` or the input list during resume. Start new runs for b3; b2 execution
contracts cannot resume after a code update.

The progress bar counts processed hours. `SUCCESS` appears only after all
chunks and the final report have been committed. Only one writer can operate
on a given output directory.

## 4. Output files and how to interpret them

```text
results/rain24_01/
  execution_contract.json
  chunk_000000/
    objects.csv
    edges.csv
    events.csv
    family_map.csv
    metadata.json
    state.json
    gap_conflicts.json
    adjacent_overlap.json
    COMMITTED.json
    identified_labels.npy       # only with --save-labels
    tracked_labels.npy          # only with --save-labels
  chunk_000006/
    ...
  summary.json
  SUCCESS
```

### Raster arrays: `.npy`

Both arrays have shape `(frames_in_chunk, y, x)`.

- `identified_labels.npy`: `int32` local object labels; restart at each frame.
- `tracked_labels.npy`: `int64` branch IDs; persist across continuation/gap links.
- Zero is background in both arrays. Missing pixels also have zero labels;
  label arrays alone do not distinguish missing data from dry/background pixels.
- Read with `np.load(path, mmap_mode="r")` for inspection without loading a
  complete saved array into RAM.

### Object catalog: `objects.csv`

One row per observed object at one frame:

```text
time,node_id,branch_id,family_id,local_label,area_cells,centroid_y,centroid_x,mean_intensity,max_intensity,precipitation_sum,y_start,y_stop,x_start,x_stop
```

`time` is the zero-based frame index across the entire run, not a Unix timestamp;
use each chunk's `metadata.json -> timestamps` for interval-ending timestamps.
For undated index-only input, timestamps are null; `sample_indices` gives the
global source sample indices instead. `time_verified` distinguishes decoded
calendar evidence from an assumed hourly timeline or an explicit date origin.
Centroids are rainfall-weighted grid coordinates. Intensities are mm/h;
`precipitation_sum` is the sum of pixel rain rates, not an area-integrated volume.
`area_cells` counts retained wet pixels; bounding-box stop indices are exclusive.
Physical area needs actual cell areas.

### Relationships: `edges.csv` and `events.csv`

```text
# edges.csv
time,parent_node_id,child_node_id,score,event,gap,raw_iou,advected_iou,distance,parent_coverage,child_coverage

# events.csv
event_id,time,event,parent_node_ids,child_node_ids
```

Each edge connects observed node IDs; `time` is the child frame. Edge events are
`continue`, `gap_continue`, `split`, `merge`, or `complex`.
`gap` counts absent object frames between endpoints; ordinary adjacent edges have zero.
`distance` is predicted-centroid error in grid cells. IoU/coverage are fractions.
`score` is a matching score, not a probability. Under overlap-ranked matching,
some accepted overlap links can be below the ordinary score threshold.
Event parent/child ID lists are semicolon-separated, e.g. `12;13`.
An event with multiple endpoints may emit several rows in `edges.csv`.

### IDs and final families

| Identifier | What it identifies | When it changes |
|---|---|---|
| `node_id` | One object observation | Every new observation |
| `branch_id` | One uninterrupted one-to-one track | New birth or split/merge/complex-event branch |
| `family_id` | Branches joined by split/merge relationships | Later merges can join previously separate families |

Continuation/gap edges preserve a branch. Split/merge ends old branches and
creates new ones. IDs are not recycled within the sequence; filtering a plot
must not compact or renumber them. Independent runs may allocate different IDs.

`family_map.csv` contains `branch_id,family_id` as known when that chunk was
written. For complete-lifecycle analysis, export final roots:

```bash
python scripts/canonicalize_stream_families.py results/rain24_01 \
  --output-dir results/rain24_01_final_families
```

**Output:** combined `objects.csv` with an added `canonical_family_id` column,
`family_map.csv` with `branch_id,canonical_family_id`, `summary.json`, and
`SUCCESS`. Requires a completed source stream and a new destination. It retains
original node/branch/family columns, does not rewrite source chunks, and does
not copy edge tables or label arrays.

### Metadata and restart files: JSON

- `metadata.json`: interval start/end times, source files, units, and grid attributes.
- `state.json`: active/dormant objects, motion history, ID counters, and family state.
- `COMMITTED.json`: chunk frame range, timings/counts, raster hashes, and payload hashes.
- `execution_contract.json`: planned inputs, preset, source hashes, and dependency versions.
- `summary.json`: completed frame/chunk counts, nodes/edges/events/gaps, and stage timings.
- `gap_conflicts.json` / `adjacent_overlap.json`: matching-decision records for optional diagnosis.

## 5. Python functions for already prepared hourly rain

These functions accept rainfall arrays, **not cumulative RAINNC**. Convert
accumulation to hourly rain rates before calling them. They do not read NetCDF
or write files automatically. Use an explicit preset to match the main runner:

```python
import json
from pathlib import Path
import numpy as np

preset = json.loads(Path("configs/rainfall_lineage_4km_1h_reference_v1.json").read_text())
rain = np.load("rain_mm_h.npy")  # Floating (time, y, x) rates in mm/h.
radius = preset["identification"]["bridge_radius_cells"]
yy, xx = np.ogrid[-radius:radius + 1, -radius:radius + 1]
structure = xx * xx + yy * yy <= radius * radius
```

### `identify()` and `identify_frame()`

**Use:** group rainfall pixels into objects for a cube or a single frame.

```python
from step import identify, identify_frame

labels = identify(
    rain, morph_structure=structure, threshold=1.0,
    min_size=1, workers=4, valid_mask=np.isfinite(rain),
)
frame_labels = identify_frame(
    rain[0], morph_structure=structure, threshold=1.0,
    min_size=1, valid_mask=np.isfinite(rain[0]),
)
```

**Returns:** `identify()` produces an `int32 (time,y,x)` array;
`identify_frame()` produces an `int32 (y,x)` array. Positive values are local
labels and zero is background. Labels are deterministic per frame; they are
not track IDs. `min_size` counts retained wet cells, not dilated area.

**Notes:** inputs must be 3D/2D respectively, with same-shaped validity masks.
NaNs are invalid pixels. Threshold units follow the input, so use mm/h
consistently. Defaults are threshold=0 and a 3x3 footprint; passing the explicit
1-mm/h/disk settings above is necessary for the reference.
`workers` parallelizes whole frames on fork-capable hosts.

### `track_with_graph()`

**Use:** track labeled objects using their masks and rain intensities.

```python
from step import track_with_graph

tracking = dict(preset["tracking"])
tracking["km"] = tracking.pop("max_displacement")
branches, graph, state = track_with_graph(
    labels, rain, sequence_id="member0_1996", return_state=True, **tracking,
)
```

**Returns:** `branches` is an `int64 (time,y,x)` branch-ID raster;
`graph` is a `TrackGraph` holding `objects`, `edges`, `events`, and
`family_map`; `state` is a `TrackingState` for the next chunk.
Without `return_state=True`, the result is `(branches, graph)`.

**Notes:** labels and rain must have identical 3D shapes. `km` is a legacy
parameter name meaning grid cells per frame, not physical kilometers.
Substantive overlap is another candidate route, so this is not a strict bound
on all accepted centroid displacements. `workers` and `phi` are accepted for
compatibility; they do not parallelize tracking or enter its current score.
Bare tracking defaults differ from the reference preset, so pass `tracking`
as above. Integer time starts at zero and resumes at `state.last_time + 1`.

### `save_tracking_state()` and `load_tracking_state()`

**Use:** persist state and continue a sequence without restarting IDs.

```python
from step import save_tracking_state, load_tracking_state

save_tracking_state(state, "after_chunk.json")
resume_state = load_tracking_state("after_chunk.json")
next_branches, next_graph, next_state = track_with_graph(
    next_labels, next_rain, state=resume_state, sequence_id="member0_1996",
    return_state=True, **tracking,
)
```

**Output:** a JSON checkpoint written by atomic rename; loading returns a
`TrackingState`. The next graph covers the next call only, not all historical
observations. Keep prior tables and the final family mapping.

**Notes:** chunks must be consecutive and use the same sequence, grid, and
tracking configuration. State does not store previous CSV/raster files or
check raw-input files. Use the main runner's `--resume` for the complete
input/output integrity workflow.

### `Tracker.update()`, `TrackingState`, and `TrackGraph`

For custom readers that deliver one frame at a time:

```python
from step import Tracker

tracker = Tracker(
    sequence_id="member0_1996", **preset["tracking"],
)
branch_frame, frame_graph = tracker.update(0, labels[0], rain[0])
```

**Returns:** an `int64 (y,x)` branch raster and a one-frame `TrackGraph`.
Continue with increasing consecutive frame indices on the same tracker;
its `tracker.state` is the resumable `TrackingState`.

`TrackingState` stores the last frame index, active/dormant terminals, next-ID
counters, grid/configuration, and family unions. `TrackGraph` exposes observed
objects (`graph.objects`, also `graph.nodes`), relationship edges, event
records, and the current branch-to-family mapping. These are in-memory
dataclasses, not files. Prefer loading state rather than manually editing counters.

**Saving Python API results:** use labels, branches, graph, and state from the
same processed chunk. Run from the source checkout and use the CSV helper to
obtain the schemas in section 4:

```python
import sys
sys.path.insert(0, "scripts")
from run_npy_validation import write_graph

out = Path("results/api_chunk_01")
out.mkdir(parents=True, exist_ok=False)
write_graph(graph, out)
np.save(out / "identified_labels.npy", labels)
np.save(out / "tracked_labels.npy", branches)
save_tracking_state(state, out / "state.json")
```

This saves that call's tables/arrays/state; it does not create a reference-runner
execution contract or `SUCCESS` marker.

## 6. Auxiliary tools and optional measurements

You do not need these for the cumulative-RAINNC workflow above.
Usage, input requirements, output formats, and notes are grouped in
**[Auxiliary functions](docs/AUXILIARY_FUNCTIONS.md)**:

- installation smoke check, regression tests, and saved-stream equivalence audits;
- plots and performance diagnostics;
- 0.5/0.1-mm/h seeded rainfall expansion and area/rain-volume measurements;
- persistent identification workers and tiled morphology;
- preprocessed NumPy input, catalog finalization, and the legacy `track()` API;
- original STEP comparisons and historical research protocols.

The recommendation remains the 1-mm/h preset. Lower-threshold envelopes are
optional measurement layers and are not enabled as matching inputs by that preset.
