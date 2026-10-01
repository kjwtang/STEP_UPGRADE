# Rain-only local validation

Use a separate output directory for each run. Record the Git commit and the
configuration. Original cumulative input files are never modified.

## Input contract

- A snapshot named `cstm_d01_YYYY-MM-DD_HH:MM:SS.rain` contains NetCDF RAINNC in
  mm, either 2-D or with a singleton Time dimension. Filename time is authoritative.
- A 72-hour window requires 73 consecutive hourly cumulative snapshots. The
  first snapshot is the predecessor; output timestamps denote interval ends.
- Preparation fully reads each file and checks finite/nonnegative data, constant
  grid attributes, no accumulation decreases, and unchanged file size/mtime.
- Latitude/longitude and internal Times are unnecessary for grid-space tracking.
  Geographic plots later require a matching coordinate reference, verified for
  dimensions, projection and crop. Its reference year is separate from rain year.
- `prepare_rain_only.py` is limited to 72 hours and holds the bounded window in
  memory. For multi-year production, use streaming ingestion instead of whole-year
  preparation. File stability is checked; source checksum comparison is separate.

## Reproduce the window

Replace DATA_DIRECTORY and choose a fresh output path:

```bash
python -u scripts/prepare_rain_only.py DATA_DIRECTORY \
  --start 1996-06-01T00:00:00 --hours 72 --output-dir results/june_input

python -u scripts/validate_real_data.py results/june_input/rain_mm_h.npy \
  --rain-kind rate --units mm/h --dt-hours 1 --hours 72 --crop 0 \
  --threshold 1 --bridge-radius 9 --workers 2 --tau 0.35 \
  --max-displacement 20 --chunk-frames 6 --stream --save-labels \
  --sequence-id june1996 --output-dir results/june_baseline
```

For the candidate, use the same command with a fresh output directory and append
`--gap-conflict-policy endpoint_overlap --adjacent-policy overlap_first
--velocity-reset off --score-policy coherent_adjacent`.

```bash
python scripts/compare_stream_tracks.py results/june_baseline results/june_candidate \
  --output-dir results/june_comparison
python scripts/replay_stream_pairs.py results/june_candidate \
  --output-dir results/june_checkpoint7 --no-focus --checkpoint-every 7
```

To audit a particular adjacent pair, replace `--no-focus` with repeatable
`--pair PARENT_NODE:CHILD_NODE`, using IDs from this run's catalog.

The prepared NPY loader reports `time_verified=false` because it does not read
the preparation sidecar. Filename chronology is documented separately in
`input_metadata.json`; this does not indicate missing timestamps in preparation.
Do not claim the NPY validator independently checked source timestamps.

## 1996 June result and next decision

Full grid 1749 by 2049, 72 hourly intervals, 16,769 identical identified objects.
Baseline: 5,098 edges, 40 gaps, two split events. Candidate: 5,317 edges, 29 gaps,
24 split and 21 merge events. Candidate adds 273 pairs, removes 54 and changes
26 event labels; these totals do not establish meteorological correctness.

Four inspected lost large-object continuations score 0.3088–0.3456 against tau
0.35. All are admitted candidates but rejected by score, so this is distinct
from a spatial gate or identifier-allocation failure. The strong raw-overlap
rule does not rescue them. Review time sequences and event geometry before
changing thresholds. Of the candidate's 29 gaps, 24 involve an endpoint under
16 cells; lowering the rain threshold globally would not specifically solve
this tracking issue.

The overlap reuse optimization reproduced all saved masks, CSV catalogs and
state JSON exactly. Restoring checkpoints every seven frames reproduced all
72 masks and graph catalogs from the original six-frame stream. Keep the
scientific policies experimental while gathering additional independent cases.
