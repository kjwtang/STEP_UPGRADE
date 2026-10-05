# RCC high-resolution validation and performance testing

Historical validation runner guide. For the current frozen beta entry point and pinned installation,
use [the group beta guide](https://github.com/kjwtang/STEP_UPGRADE/blob/v0.4.0b2/docs/GROUP_BETA_0_4_0b2.md). This older runner is not the reference resume interface.
Use a separate upgraded environment, not the original STEP environment. Commands run from the
repository root; install dependencies on a login node and compute on a compute node:

```bash
python -m pip install -r requirements-validation.txt
```

## RCC CIMP_1hr: original hourly CSTM directory

The inspected data has one file per hour, cumulative `RAINNC(south_north,west_east)` in mm,
and a 1010x1634 full domain. This dedicated reader loads cropped rainfall and coordinates,
not 3D meteorological fields. It recursively finds `cstm_d01_*.nc`, sorts filename timestamps,
and checks each against internal Times. `--start 0 --hours 24` reads the earliest 25 snapshots
to produce 24 hourly increments; the first snapshot is only the difference baseline.
Each chunk includes its preceding snapshot.

Historically, update the development checkout with `git pull --ff-only` on the login node,
then run on the compute node. Do not update a pinned beta checkout during an active job.

```bash
python scripts/validate_real_data.py \
  /project2/moyer/kjwtang/tmp/CIMP_1hr/2005.07 \
  --wrf-cumulative --inspect

/usr/bin/time -v python -u scripts/validate_real_data.py \
  /project2/moyer/kjwtang/tmp/CIMP_1hr/2005.07 \
  --wrf-cumulative --rain-kind cumulative \
  --hours 24 --crop 600 --workers 4 --chunk-frames 6 \
  --sequence-id present_2005_rainnc_validation \
  --stream --output-dir results/benchmark_rainnc_600
```

After inspection, use crop=0 and `results/benchmark_rainnc_full` for the full native grid.
For small-sample whole/chunk equivalence, remove `--stream`, use `--crop 300`, and select a new output directory.

Missing hours, duplicate times, filename/Times mismatch, grid changes, and cumulative decreases
stop processing. Negative-difference tolerance defaults to zero. Set `--negative-tolerance-mm`
only with evidence of floating-point rounding; clipped small-negative pixel counts are recorded.
Restarts/cumulative-bucket rollback are not automatically converted to zero rain.
Intervals not recoverable from adjacent accumulations require original data or restart records.
A mismatched `wrfout_ref` year triggers a warning and retains provenance in chunk metadata;
even matching Times/filenames do not prove the source was not copied or mislabeled.
This validates RAINNC only, not an unavailable RAINC convective component.

## 1. Inspect inputs first

```bash
python scripts/validate_real_data.py /PATH/TO/rainfall.nc --inspect
```

The generic reader accepts a single NetCDF file or `(time,y,x)` NPY, not a directory;
the dedicated WRF reader above handles directories. Verify variable, dimensions, time, and units.
`PRECIP` and `time y x` below are placeholders to replace with inspected names.
Raw cumulative RAINNC/RAINC requires differencing and restart handling before generic use.
Already processed hourly PRECIP(mm) uses `--rain-kind interval`, without another difference pass.
Use `--units mm` only when units are known but metadata is absent. NPY also requires
`--dt-hours 1 --units mm`; chronological continuity cannot be automatically verified.

## 2. Small-sample consistency validation

```bash
python -u scripts/validate_real_data.py /PATH/TO/rainfall.nc \
  --variable PRECIP --dims time y x --rain-kind interval \
  --start 0 --hours 72 --crop 300 --workers 4 \
  --chunk-frames 24 --sequence-id present_2005_validation \
  --output-dir results/check_300_72h --gif
```

`--hours` counts frames and equals hours only for hourly data. Historical defaults
threshold=1 mm/h, radius=9 cells, displacement=20 cells/frame are calibration starting points.
Adjust grid parameters to maintain physical scales when resolution changes. These commands
expand native-resolution coverage; they do not downsample.

Require all `summary.json` checks=true: whole/chunk identification and labels, final family mappings,
nodes/edges/events, unique nodes, per-frame unique branches, valid forward-time edge references.
Plot numbers are branch IDs; split/merge creates new branches intentionally. Colors repeat after 256.
`statistics.png` includes counts, area, observed branch spans, and within-branch speeds.
Spans include gaps and are truncated by time/crop boundaries and split/merge; they are not true storm lifetimes.
`frame_*.png` includes sampled frames, first event-type neighborhoods, and first gap endpoints.
GIFs are optional and may be slow/memory-intensive. Coordinates are cropped grid cells, not latitude/longitude.

## 3. High-resolution streaming benchmark

Begin with 24 frames on a native 600x600 crop, reading six frames per chunk. Test 1000x1000
in another process, then full domain with `--crop 0`. Use a distinct output directory each time;
existing directories are not overwritten.

```bash
/usr/bin/time -v python -u scripts/validate_real_data.py /PATH/TO/rainfall.nc \
  --variable PRECIP --dims time y x --rain-kind interval \
  --start 0 --hours 24 --crop 600 --workers 4 --chunk-frames 6 \
  --sequence-id present_2005_benchmark --stream \
  --output-dir results/bench_600
```

Streaming reads, identifies, tracks, saves, and reloads checkpoints per chunk without retaining
all rainfall rasters. Boundary times are checked; frames are not silently skipped. Tracking within
a sequence is serial, while identification is parallel. Outputs include chunk tables/checkpoints;
`--save-labels` additionally saves identification/tracking arrays. Previews cover the first chunk;
benchmarks usually omit GIFs. Size resources using smaller-domain memory and object counts first.
`SUCCESS` means the runner completed. Preserve failed outputs for diagnosis; this historical runner
does not automatically restart in the same output directory.

Outputs:

- `performance.csv`: per-chunk read/ID/tracking/output/plot times, object counts, memory peaks.
- `performance.png`: total stage timings and memory trends.
- `summary.json`: core time per frame, million grid cells per second, object/event/gap counts.
- `frame_statistics.csv`, `statistics.png`: object counts, mean/max rainfall, wet/missing fractions.
- `chunk_*/`: node/edge/event/family tables, time/unit metadata, checkpoints.

Process-tree RSS is sampled every 0.2 seconds; fork-shared pages can be counted repeatedly,
so it is approximate, not exclusive physical memory. `parent_peak_rss_mb` excludes ID workers.
Cross-check `/usr/bin/time -v` and post-job
`sacct -j JOBID --format=JobID,Elapsed,TotalCPU,MaxRSS,State`.
Parent CPU time also excludes ID workers. Runner wall time includes output and first-chunk plots,
but excludes dependency import and final performance-summary plots; outer time covers the full command.
Core performance uses identification+tracking separately. Restricted memory queries populate
`memory_sampling_errors`; unavailable values are null, and partial samples are not a complete tree peak.

Streaming does not automatically establish whole-run equivalence or scientific validity;
small-sample mode checks consistency. Object extraction, pixel overlap, and global assignment
can slow sharply with object count. Chunked input does not guarantee arbitrary domains fit 32 GB.
Inspect dense-rainfall cases; do not extrapolate an entire JJA solely from grid-cell count.
This historical runner needs external planning for cross-month production. The generic reader accepts
one multi-time file; only the dedicated WRF option accepts original hourly CSTM directories.
Do not mix climate years/members in one input directory. The beta reference runner provides its own continuous plan.

Record commit, interval, threshold, radius, displacement, workers, and chunk size when scaling up.
Where needed, compare workers=1/4 and chunks=3/6 with everything else fixed.
Consistency checks do not replace physical inspection of false links, breaks, split/merge, and gaps.
