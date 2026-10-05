# Multi-time NetCDF input validation: 2026-10-05

Target release: v0.4.0b3. The identification/tracking source and the
`rainfall_lineage_4km_1h_reference_v1` preset are unchanged from v0.4.0b2.
This validation concerns ingestion and execution equivalence, not a new
meteorological tracking accuracy assessment.

## Automated checks

- Full suite: 206 passed, one skipped (original-STEP checkout not provided).
- New reader regressions: 30 passed. Coverage includes cumulative/interval/rate
  classification, conversion to mm/h, ambiguous fields, numeric hourly defaults,
  decoded/string/WRF character times, multi-file ordering, boundary duplicates,
  missing hours, resets, missing-pixel barriers, grid/coordinate changes, bounded
  reads, changed-input rejection, and checkpoint resume.
- Installation smoke: all ten checks true; 12 frames, 11 nodes, 10 edges and one
  gap edge. Snapshot, single-file and multi-file runs match exactly, including
  changed chunk sizes and workers after a controlled pause.
- `pip check`: no broken requirements. Existing NumPy/netCDF runtime and xarray
  deprecation warnings remain in this pinned environment; they do not change
  the recorded equality results. Run the smoke check on the target RCC node.

Environment: Python 3.13.5, macOS ARM64; dependency versions from
`requirements-beta.txt` (unchanged top-level pins from b2).

## Local real-data equivalence

Input: 13 unmodified cumulative RAINNC snapshots from
1996-01-01 00:00 through 12:00, hourly, full 1749 x 2049 grid, DX=DY=4000 m.
The local snapshots were repackaged into one multi-time file and two multi-time
files solely for comparison. Original data was read-only; datasets/results are
not committed to GitHub.

| Execution | Frames | Nodes | Edges | Split events |
|---|---:|---:|---:|---:|
| Original hourly snapshot directory, 6-hour chunks | 12 | 2034 | 568 | 3 |
| One multi-time NetCDF, 4-hour chunks | 12 | 2034 | 568 | 3 |
| Two multi-time files, pause after 4 hours, resume in 3-hour chunks | 12 | 2034 | 568 | 3 |

Both repackaged executions passed all eight complete comparisons against the
snapshot control: canonical objects, full edges/scores, events, final family
roots, saved timestamps, grid/units/cadence, identification raster hashes and
tracked raster hashes. The second file's first sample supplied an ordinary
cross-file cumulative difference; it was not dropped or treated as a new sequence.

## Interpretation and input boundaries

Calendar times verify hourly continuity across files. Local numeric indices
cannot independently reveal a missing monthly file; their declared/default
one-hour cadence and ordered-file continuity are recorded as assumptions.
The reader never classifies precipitation from its numerical monotonicity.
Cumulative resets are rejected, not clipped; an hourly predecessor is reused
across processing chunks and tracking state is restored from committed JSON.

Only requested rainfall slices are loaded. The time/shape/file metadata plan is
established before execution. Resume validates the same planned input extent,
source code, environment and output contract. Changing the planned file list,
appending hours, or updating a paused checkout is not a supported resume action.
Preserve the original checkout to resume older beta runs.
