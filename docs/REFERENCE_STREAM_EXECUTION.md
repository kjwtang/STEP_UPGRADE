# Restartable execution of the frozen rainfall reference

This is an **execution-layer extension**, not a new storm definition. The
ranked-v2 tracking core and 1-mm/h / radius-9 identification reference remain
unchanged. Defaults of `identify` and `track_with_graph` are unchanged.

## Execution contract

`scripts/run_reference_stream.py` loads the named reference JSON directly.
Its specialized reader requires full-grid, hourly filename timestamps,
cumulative `RAINNC` in millimetres, and constant WRF projection attributes with
DX=DY=4000 metres. N intervals need N+1 snapshots, including the predecessor.
Missing observations, accumulation decreases and changed grids fail; there is
no automatic reset repair or extrapolation of the first interval.

The runner scans the directory once, processes bounded batches, parallelizes
independent-frame identification, and carries one ordered tracker state forward.
It does **not** independently track spatial tiles or months and then renumber
their tracks. Absolute tracker time is relative to the beginning of this planned
sequence. Metadata preserves each interval's actual filename timestamp.

The execution contract records the complete selected snapshot plan (path,
timestamp, byte size and nanosecond mtime), full preset, grid, sequence namespace,
relevant source SHA-256 digests and dependency versions. These guards supplement
the tracker's own parameter/shape checks, including the upstream ID threshold.
File size/mtime is a transfer/change detector, **not a cryptographic digest of
input content**. Filename time and projection attributes are not independent
proof of geographical placement, source-year provenance, or rain accuracy.

Resume requires the same plan, namespace, scientific settings, source code,
dependencies and label-saving choice. Workers, batch length and the verified
fresh/persistent executor can change at a committed boundary. A longer plan or
a changed data location is not silently accepted as a continuation. Keep the
checkout unchanged during a run and its restarts; use a pinned worktree.

## Publication and interruption behavior

Each batch writes graph catalogs, input metadata, diagnostic decisions and the
tracking checkpoint into a hidden staging directory. Optional label arrays are
diagnostics, disabled by default. `COMMITTED.json` records file SHA-256 digests,
per-frame label/raster hashes, counts and timing, then the whole directory is
published as `chunk_000000`, etc. Graph and checkpoint cannot be published as
separate completed chunks.

Resume verifies every committed payload and contiguous frame extent, then loads
only the most recent checkpoint. Hidden `.pending` directories are ignored and
preserved as interruption evidence; they are never treated as processed frames
or automatically deleted. A POSIX advisory lock rejects a second writer to the
same output and is released if its owning process exits. The sequence is only
successful after every requested frame is committed and the final contract is
rechecked. `PAUSED.json` is historical metadata; on successful resumption it is
marked resolved. `SUCCESS` is the completion authority, not a 100% tracking bar.
Resumed ETA uses only frames processed in the current invocation as its speed
sample; the percentage still includes previously committed frames. Initial ETA
is unavailable until the first new frame completes.

These tests cover process interruption/publication failure. Atomic directory
rename is not a claim of proven power-loss durability on every filesystem, and
the advisory lock must be supported by the deployment filesystem. Verify this
on RCC before production. Output records are integrity checks, not authenticated
signatures against deliberate tampering.

## Run and restart

Install the package in the active validation environment first; scripts do not
guess another environment or silently import a different STEP checkout.
The specialized POSIX runner uses Python >=3.8 APIs; the group beta's installation
floor is >=3.12 with its pinned dependencies, and its reference environment is
Python 3.13.5. See the [group beta guide](GROUP_BETA_0_4_0b1_ZH.md) for installation.

```bash
python -m pip install -e . -r requirements-validation.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

python -u scripts/run_reference_stream.py DATA_DIRECTORY \
  --start PREDECESSOR_SNAPSHOT_INDEX --hours 2208 \
  --sequence-id historical_member0_1996_jja_reference_v1 \
  --workers 8 --chunk-frames 24 --id-executor fresh \
  --output-dir results/reference_jja

# Same command, same plan and output, with --resume after an interruption.
python -u scripts/run_reference_stream.py DATA_DIRECTORY \
  --start PREDECESSOR_SNAPSHOT_INDEX --hours 2208 \
  --sequence-id historical_member0_1996_jja_reference_v1 \
  --workers 8 --chunk-frames 24 --id-executor fresh \
  --output-dir results/reference_jja --resume
```

`--stop-after-chunks N` is a reproducible controlled pause, not completion. It
counts new batches in the current invocation. Do not confuse the predecessor
index with a rain interval's ending timestamp. A June-1 00:00 predecessor yields
the first interval ending June-1 01:00. A 2208-hour run reaches September-1 00:00,
which is the end of the last August-31 hourly interval.

## Optional persistent identification workers

`step.parallel_identification.FrameIdentifier` reuses one forked pool and bounded
shared wet/valid masks plus disjoint output slots. Workers receive frame indices
and return only completion indices, rather than pickled full label arrays.
Every batch completes before its buffers are overwritten; returned labels are
independent copies. Existing morphology, missing-data barriers, size filtering,
global connectivity and deterministic relabeling are reused exactly.

Use `--id-executor persistent` to test it. Create the pool before application
threads; this runner has no asynchronous NetCDF reader thread and never forks a
new ID pool while a persistent pool is active. Pool closure is explicit on
normal completion, pause and exception. This is not a parallel tracker or an
implemented I/O-prefetch pipeline. Shared buffer capacity increases with batch
length; shared pages must not be counted as separate physical RAM per worker.

Both executors remain available. Do not assume persistent is faster merely
because fewer pools are created. Benchmark prepared input forward and backward,
including startup/shutdown and enough batches to exercise reuse:

```bash
python -u scripts/benchmark_id_batches.py PREPARED_MM_H.npy \
  --frames 72 --batch-frames 6 24 --workers 4 8 --repeat 2 \
  --output-dir results/id_batches_forward
python -u scripts/benchmark_id_batches.py PREPARED_MM_H.npy \
  --frames 72 --batch-frames 6 24 --workers 4 8 --repeat 2 --reverse \
  --output-dir results/id_batches_reverse
```

This measures warm prepared identification only; it excludes ingestion, tracking,
hash computation and output. An end-to-end run is a separate comparison. Laptop
tests do not establish RCC 16/32-worker scaling or physical job-memory limits.

The completed 72-frame full-grid prepared-input benchmark ran forward and
reverse, with two repeats each. Median identification seconds across four
measurements per configuration, including pool startup/shutdown:

| Batch frames | Workers | Fresh pool | Persistent pool |
|---:|---:|---:|---:|
| 6 | 4 | 12.573 | 12.221 |
| 6 | 8 | 7.021 | 6.530 |
| 24 | 4 | 9.710 | 9.590 |
| 24 | 8 | 5.197 | 5.046 |

Every identified frame matched the serial reference hash. Persistent reuse
improved this stage by approximately 1–7%, not a demonstrated substantial
end-to-end speedup. A separate single-batch 24-frame experiment was slower
with persistent workers. The default therefore remains `fresh`; the optional
executor is not a scientific change or an automatic performance promotion.

## Evidence and limitations

The initial month experiment processed **744 full-grid 1996 July intervals**
(1749 x 2049), pausing after 48 frames and resuming with changed workers and batch
size. Against the prior validation runner, all 175,329 object measurements and
immutable IDs, 53,732 full edges/scores, 392 event records, and final family roots
match. Its first 168 frames match the previously saved identification and
tracking rasters exactly. The old monthly control did not save full-month
rasters, so full-month pixel equality was unavailable in that initial comparison.

The completed **2208-hour 1996 JJA** test contains 504,350 objects, 154,046 edges,
352,793 branches, and 1,238 event records (688 split / 542 merge / 8 complex).
A fresh-worker 24-frame control was compared with persistent workers paused
after 720 frames at the June/July handoff, then resumed in 32-frame batches.
Every frame's ID and tracked-raster hash, every immutable object field, full
edge including score/overlap fields, event, saved timestamp/grid/unit contract,
and final-known family root match. This is exact execution equivalence on one
year's summer, not new algorithm calibration or multi-year scientific truth.

A real SIGINT test interrupted a 72-frame stream after its first six-frame
chunk published. Resume with changed workers/executor/batch length reproduces
the prior frozen June-1 reference exactly, including all raster hashes and
catalogs. A synthetic publication failure also leaves an ignored/preserved
hidden staging directory and reprocesses its unpublished frames.

The JJA final checkpoint is 8,449,545 bytes and still retains 352,793 historical
branch parents. Checkpoint history is **not bounded** by batch size. All saved
checkpoints total 429,975,597 bytes in the 24-frame control versus 341,117,315
bytes with the changed boundary schedule. Different totals reflect different
snapshot frequency and object activity; they are not a ten-year space forecast.
No worker-summed physical-memory bound is established by these laptop tests.

An additional untouched winter window covers 120 intervals ending February-27
01:00 through March-3 00:00, including February-29. A seven-frame fresh-worker
control matches a 13-frame persistent-worker run paused after 39 frames and
resumed with five workers instead of eight. All raster hashes and full catalogs
match: 22,806 objects, 7,444 edges, 20 gap edges, 53 event records. JJA plus this
window comprise 2,328 distinct intervals from the same 1996 year; repeated
execution is not additional independent validation data.

The final synthetic/diagnostic suite has 138 passed and one skipped test; the
skipped test requires a separately pinned original STEP checkout, unavailable
locally. `pip check` reports no broken declared requirements. A native NumPy /
NetCDF4 binary-size warning and NumPy-2.5 deprecation warnings remain in this
environment; passing local tests does not resolve binary-ABI compatibility or
certify the RCC environment. Test a clean pinned deployment environment there.

The new audit tool normalizes family IDs to the last checkpoint before catalog
comparison; per-chunk as-of roots can legitimately differ with batch boundaries.
It checks node uniqueness, branch non-reuse and connectivity, forward/gap time,
event references, next-ID counters and final family roots. Exact raster hashes
are compared when available; null means unavailable, not passed:

```bash
python scripts/audit_reference_stream.py results/CANDIDATE \
  --reference results/CONTROL --output-dir results/stream_equivalence
```

No observational ground truth, new climate-member validation, full ten-year
production qualification or 0.1-mm/h envelope validation is established here.
Historical union-find roots still grow with the sequence and are serialized
into every checkpoint. Output hashing and retained checkpoint archives cost
time and space; retention/pruning requires a separately verified design.
Final-known families can be exported with `canonicalize_stream_families.py`
without rewriting immutable branch IDs or source catalogs.
