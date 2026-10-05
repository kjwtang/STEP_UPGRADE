# Auxiliary functions and optional workflows

These tools supplement the main cumulative-RAINNC workflow. Run CLI commands
from the installed beta checkout; use a new output directory for each new run.
Examples use hourly rain rates in mm/h where an array, rather than a cumulative
file directory, is requested. Python examples use prepared arrays and the
`preset`/`structure` variables introduced in README section 5.

## 1. Installation and result checks

### `scripts/beta_smoke.py`

**Use:** test the installed beta without downloading research data.

```bash
python -u scripts/beta_smoke.py --output-dir results/beta_smoke_01 --workers 2
```

**Input:** no external data; `--workers` must be positive. The installed package
version must match the script's beta version.

**Output:** generated NetCDF snapshots under `input/`, whole-run `control/`,
paused/resumed `resumed/`, `summary.json`, and `SUCCESS`.
Expected summary: 12 frames, 11 nodes, 10 edges, one gap edge, eight exact checks.
The test crosses a month boundary and changes chunk size/worker count on resume.
Existing output directories are rejected.

### Regression tests

```bash
python -m pytest -q
```

**Output:** terminal pass/fail/skip summary; pytest's local cache.
Tests cover identification, graph/ID rules, split/merge/gaps, saved state,
worker/chunk equivalence, input validation, envelopes, and the English-only text policy.
The original-STEP comparison test requires a separately pinned original checkout
and can be skipped when that reference is absent.

### `scripts/audit_reference_stream.py`

**Use:** check saved graph/ID structure, or compare two executions of the same
input/preset with different workers/chunks.

```bash
python scripts/audit_reference_stream.py results/rain24_01 \
  --output-dir results/rain24_structure

python scripts/audit_reference_stream.py results/rain24_other_workers \
  --reference results/rain24_01 \
  --output-dir results/rain24_equivalence
```

**Input:** completed stream directories with `SUCCESS`, chronological chunk
tables, metadata, and checkpoints. The reference must describe the same input
interval and grid. Reference-runner chunks carry hashes even without saved arrays.

**Output:** `summary.json`, `REPORT.md`, and `SUCCESS`.
Checks include unique/monotonic IDs, forward-time edges, branch continuity,
event endpoints, final family roots, full scores/tables, metadata, and raster hashes.
A null comparison means evidence was unavailable, not that the check passed.
`--prefix-reference SHORT_RUN` compares an earlier prefix's timestamps and
raster hashes; it does not establish whole-graph equality for different time spans.

## 2. Plots and performance diagnostics

### `scripts/validate_real_data.py`

**Use:** run a short diagnostic case with plots, statistics, and whole/chunk
comparison, or run a streaming timing/memory benchmark. This is a separate
diagnostic runner: it reads and reruns the case, not a plot-only command for
reference-runner output.

For filename-timed RAINNC-only files, specify both reader modes. The following
tracking settings reproduce the reference policies explicitly:

```bash
python -u scripts/validate_real_data.py /PATH/TO/ONE_SEQUENCE \
  --wrf-cumulative --rain-kind cumulative \
  --wrf-time-source filename --wrf-grid-source attributes \
  --hours 12 --crop 300 --workers 4 --chunk-frames 6 \
  --threshold 1 --bridge-radius 9 --min-size 1 \
  --tau 0.35 --max-displacement 20 --max-gap 1 \
  --gap-tau 0.45 --gap-ambiguity 0.05 --event-overlap 0.1 \
  --gap-conflict-policy endpoint_overlap --gap-min-pixels 16 \
  --adjacent-policy overlap_ranked --score-policy coherent_adjacent \
  --velocity-reset off --family-map-scope observed \
  --sequence-id member0_diagnostic \
  --output-dir results/diagnostic_12h
```

**Input:** cumulative WRF directory with the options above, or one prepared
NetCDF/NPY cube with its variable/dimensions/units configured. For prepared NPY
rates add `--rain-kind rate --units mm/h --dt-hours 1` instead of WRF options.
`--inspect` prints input metadata without running ID/tracking.

**Output, non-stream mode:** object/edge/event/family CSV tables, label NPY arrays,
JSON settings/metadata/checkpoints/summary, `statistics.png`, sampled
`frame_*.png`, and `SUCCESS` after checks. `--gif` adds an animated preview.
The summary records whole/chunk equivalence and structural checks.

**Benchmark mode:** add `--stream`; use `--crop 0` for the full domain.
Outputs include `performance.csv`, `performance.png`,
`frame_statistics.csv`, `statistics.png`, chunk tables/checkpoints, summary,
and first-chunk previews. Add `--save-labels` to retain per-chunk NPY arrays.
Streaming benchmarks do not run a separate whole-sequence equivalence control
or provide the reference runner's transactional `--resume` workflow.

**Notes:** `--hours` counts frames; cropped coordinates are grid indices.
`--grid-km 4` enables nominal uniform-grid area/speed statistics when that
assumption is appropriate. Process-tree RSS is sampled and can count shared
pages more than once. Use job accounting for allocation decisions. Full plotting,
GIFs, and saved labels increase I/O and memory use.

## 3. Optional lower-threshold rainfall measurements

The examples below expand a saved seed frame without changing its identification
or tracking graph. They are not automatically run by the reference runner.

### `expand_seeded_envelope()`

**Use:** assign connected weak rain to existing seed labels.

```python
from step.envelope import expand_seeded_envelope

envelopes = expand_seeded_envelope(
    rain_frame, seed_labels, envelope_threshold=0.1,
    mode="partition", connectivity=8,
)
```

**Input:** a 2D rain-rate array and same-shaped nonnegative integer seed labels.
Every seed pixel must be finite and >= the envelope threshold.
Use local seed labels within the int32 range; zero is background.

**Returns:** an `int32 (y,x)` label array. `partition` preserves seed label
values and assigns weak rain by marker watershed. Components without a seed
remain zero. Original seed arrays are not modified.

**Notes:** `mode="system"` instead returns weak-rain connected components
containing at least one seed; those component labels are not original storm IDs.
Neither mode crosses dry/missing pixels. Connectivity is 4 or 8.
An optional `max_distance_km` requires `spacing_km=(dy,dx)` on a uniform grid.
It caps Euclidean distance to any seed before partitioning, not own-label
distance or distance along a wet path. Partition mode needs scikit-image
(included in the beta requirements).

### `identify_envelope()`

**Use:** create seed labels and their envelopes together for a 2D frame.

```python
from step.envelope import identify_envelope

cores, envelopes = identify_envelope(
    rain_frame, core_threshold=1.0, envelope_threshold=0.1,
    mode="partition", min_core_cells=1, connectivity=8,
)
```

**Returns:** `(cores, envelopes)`, both same-shaped `int32` arrays.

**Notes:** require `0 < envelope_threshold <= core_threshold`. By default
seeds are connected threshold patches, not STEP's radius-9 grouped objects.
Pass the reference morphology footprint as `morph_structure` to use that
grouping rule, or use `expand_seeded_envelope()` for already identified seeds.
`min_core_area_km2` requires explicit positive `cell_area_km2` values.
Changing seed/minimum-size settings changes which objects enter this optional
measurement, so it is not a substitute for the full reference catalog.

### `envelope_statistics()`

**Use:** measure the seed/envelope relationship and rain contributed by the
expanded area.

```python
from step.envelope_statistics import envelope_statistics

rows, memberships = envelope_statistics(
    rain_frame, seed_labels, envelopes,
    cell_area_km2=cell_areas,
)
```

**Input:** same-shaped 2D rate/seed/envelope arrays. Every retained seed pixel
must lie in an envelope. `cell_areas` is an explicit positive scalar or 2D
array in km2; omit it if physical areas are not available.

**Returns:** two lists of dictionaries, suitable for CSV/JSON serialization.

- `rows`: `envelope_label, core_count, envelope_cells, core_cells,
  envelope_to_core_cell_ratio, peak_mm_h, core_peak_mm_h, mean_mm_h,
  rain_rate_sum, core_rain_rate_sum, noncore_rain_rate_fraction,
  touches_domain_boundary, touches_missing, envelope_area_km2,
  core_area_km2, rain_volume_rate_m3_h`.
- `memberships`: `envelope_label, core_label, core_cells`, retaining
  many-to-many associations when a grouped seed appears in several system envelopes.

**Notes:** without cell areas, physical-area and volume-rate fields are null.
Volume rate is `mm/h * km2 * 1000 = m3/h`, not accumulated storm volume.
Integrate over time separately if cumulative volume is needed. The function
does not save files or alter tracks.

## 4. Optional identification execution modes

### Persistent workers: `FrameIdentifier`

The main runner supports `--id-executor persistent` instead of `fresh`.
It reuses workers/shared masks between chunks; CSV, NPY, and checkpoint schemas
are unchanged. Compare exact results and timing on your target machine first.

For a custom prepared-array reader:

```python
from step.parallel_identification import FrameIdentifier

with FrameIdentifier(
    grid_shape=rain.shape[1:], morph_structure=structure,
    capacity=6, workers=4, threshold=1.0, min_size=1,
) as identifier:
    labels = identifier.identify(rain[:6])
```

**Input:** each batch is a `(frames,y,x)` rate cube, with 1..capacity frames
and the configured fixed grid. Optional validity masks must match the batch.

**Returns:** independent-copy `int32` label arrays; a later batch cannot mutate
previous outputs. The context manager closes the pool. Outside a `with` block,
call `identifier.close()`; after closing, further identification is rejected.
Multiple workers require fork support.

### Tiled morphology through `identify()`

```python
labels = identify(
    rain, structure, workers=4, threshold=1.0, min_size=1,
    parallel_mode="tiles",
)
```

**Returns:** the same `int32 (time,y,x)` identification format.
Tiles parallelize dilation with halos; connectivity and labeling are resolved
globally. They do not independently track spatial tiles or stitch track IDs.
Use `parallel_mode="frames"` for the usual whole-frame execution.
More workers are not an automatic speedup; choose using a local benchmark.

## 5. Prepared-array and legacy compatibility tools

### `scripts/prepare_rain_only.py`

**Use:** make a small prepared rate cube from a flat directory of cumulative
`.rain` NetCDF files when an array-based diagnostic is needed.

```bash
python scripts/prepare_rain_only.py /PATH/TO/ONE_MONTH \
  --start 1996-07-01T00:00:00 --hours 12 \
  --output-dir results/prepared_12h
```

**Output:** `rain_mm_h.npy`, `float32 (time,y,x)`, and
`input_metadata.json` with start/end times, source files, and grid metadata.

**Notes:** accepts 1..72 intervals, exact `cstm_d01_*.rain` filenames, and
a timestamp for the preceding cumulative snapshot. Unlike the main runner,
it does not scan nested folders or accept `.nc` aliases. It rejects nonfinite
inputs and accumulation decreases. It holds this short cube in memory;
it is not needed for normal cumulative-RAINNC execution.

### `scripts/run_npy_validation.py` and `scripts/finalize_catalog.py`

These are older array-based table/export interfaces. Prefer the main runner
for the frozen preset and full restart checks.

```bash
python scripts/run_npy_validation.py results/prepared_12h/rain_mm_h.npy \
  --start 0 --stop 12 --threshold 1 --bridge-radius 9 --min-size 1 \
  --workers 4 --tau 0.35 --max-displacement 20 \
  --sequence-id member0_array_demo \
  --state-output results/array_demo/state.json \
  --output-dir results/array_demo
```

**Input:** a prepared `(time,y,x)` NPY rate cube, not cumulative snapshots.
`--start/--stop` slice the input cube; `--absolute-start` sets the graph's
frame index. To continue, provide `--state-input` and identical tracking
parameters/sequence. This CLI exposes historical score-policy defaults,
not the complete reference-preset policy selection.

**Output:** identified/tracked NPY arrays, object/edge/event/family CSV tables,
`run_settings.json`, and an optional JSON state file inside the output directory.

To combine chronological parts from this older runner:

```bash
python scripts/finalize_catalog.py results/array_june results/array_july \
  --output-dir results/array_catalog
```

**Output:** combined object/edge/event/family CSV tables and
`catalog_manifest.json`; label rasters stay in their original partitions.
Inputs must carry state continuously and supply all branches/endpoints;
the final part must have the complete family map. This tool updates output
`family_id` using that map. For reference-runner streams with observed-only
maps, use `canonicalize_stream_families.py` instead, as shown in the README.

### `track()`: raster-only compatibility API

```python
from step import track

options = dict(preset["tracking"])
options["km"] = options.pop("max_displacement")
branches = track(labels, rain, phi=None, **options)
```

**Returns:** only the `int64 (time,y,x)` branch raster, not a graph or checkpoint.
The example explicitly supplies reference policies from the preset;
do not include `return_state=True`. `phi` is ignored, `workers` does not
parallelize tracking, and `km` means cells per frame.
Use `track_with_graph()` for graph/state access and continuation.

## 6. Historical comparisons and research protocols

These detailed protocols include their commands, inputs, and output tables.
They are optional; the main workflow does not require rerunning them.

- [Original STEP comparison](ORIGINAL_COMPARISON.md): pinned upstream source,
  matched rainfall/object controls, plots, stage timing/memory, and failure logs.
- [Real-data diagnostics](REAL_DATA_VALIDATION.md): additional NetCDF/NPY plotting and resource checks.
- [Reference execution](https://github.com/kjwtang/STEP_UPGRADE/blob/v0.4.0b2/docs/REFERENCE_STREAM_EXECUTION.md):
  restart/worker audits and existing long-run controls.
- [Frozen tracking configuration](https://github.com/kjwtang/STEP_UPGRADE/blob/v0.4.0b2/docs/TRACKING_REFERENCE_FREEZE.md):
  exact policies and implementation reference.
- [Seeded-envelope experiments](https://github.com/kjwtang/STEP_UPGRADE/blob/v0.4.0b2/docs/CORE_ENVELOPE_EXPERIMENTS.md):
  threshold/ownership alternatives and comparison outputs.
- [Rainfall-threshold sensitivity](https://github.com/kjwtang/STEP_UPGRADE/blob/v0.4.0b2/docs/RAIN_THRESHOLD_SENSITIVITY.md):
  changed identification versus unchanged-graph measurements.
- [Tracking-geometry sensitivity](https://github.com/kjwtang/STEP_UPGRADE/blob/v0.4.0b2/docs/TRACKING_GEOMETRY_SENSITIVITY_2026-10-05.md):
  expanded-feature tracking controls and saved differences.
