# Local candidate assessment — 1996 hourly rainfall

Tracking core tested at `0af014c`, balanced-overlap **v2**, full 1749 x 2049 grid.
The scientific default remains unchanged; this is an opt-in candidate, not a
validated convective-storm definition. See [fixed policy and acceptance gates](TRACKING_FREEZE_PROTOCOL.md).

## Real-input windows

RAINNC-only files were read without modification, differenced hourly, and checked
using filename chronology and constant grid attributes. Windows start with a
midnight cumulative snapshot; the first derived interval ends at 01:00.

| Starting snapshot | Intervals | Objects | Baseline edges | Candidate edges | Candidate gaps | Events: split / merge / complex |
|---|---:|---:|---:|---:|---:|---|
| June 1 | 72 | 16,769 | 5,098 | 5,272 | 6 | 25 / 23 / 1 |
| June 30 | 72 | 17,517 | 4,533 | 4,714 | 15 | 21 / 18 / 0 |
| July 1 | 168 | 39,862 | 11,440 | 11,809 | 41 | 40 / 36 / 1 |
| August 15 | 72 | 16,218 | 4,463 | 4,632 | 6 | 23 / 17 / 1 |

384 processed intervals include 48 overlapping intervals in June-30/July-1
windows: **336 distinct hours**, not 384 independent samples. All comparisons
preserve identified masks and measurements exactly. Structural checks pass:
forward-only graph edges, unique nodes, continuous non-reused branch IDs, and
gap lengths consistent with timestamps. The June-1 72-frame stream reproduces
all rasters and CSV catalogs when restored every seven frames instead of six.

Largest added/removed pairs and split/merge mask examples were inspected locally;
they show useful recovery of organized rain-band connections, but also remaining
moderate-overlap interruptions. For example, a removed baseline continuation in
the June-30 window has raw endpoint coverages 0.259/0.377, below the candidate's
symmetric rule. Baseline is not truth; this is a residual ambiguity, not evidence
that every rejected connection is erroneous. Some removed gaps are replaced by
adjacent lineage, but many are not: do not equate fewer gaps with improvement.

In the four candidate windows, respectively 371, 353, 802 and 277 links have zero
raw **and** translated overlap; none has both endpoints at least 100 cells. The
score rule can connect nearby small objects without mask overlap. Neither the
size cutoff nor this audit establishes physical truth. Tiny objects are retained
in identification; the 16-cell floor affects gap endpoints only.

## Worker test

24 full-domain prepared frames, same rainfall and candidate settings. Seven
configurations reproduce identical masks, tracked rasters, object measurements,
edge scores, event lists, branch IDs and family maps within each test round.
A reverse-order round with math-library thread limits gives:

| Mode | Workers | ID seconds | Sequential tracking seconds |
|---|---:|---:|---:|
| Frames | 1 | 11.57 | 3.97 |
| Frames | 2 | 6.11 | 4.07 |
| Frames | 4 | 3.23 | 4.15 |
| Frames | 8 | 1.79 | 4.04 |
| Halo tiles | 2 | 6.68 | 4.07 |
| Halo tiles | 4 | 4.13 | 4.03 |
| Halo tiles | 8 | 2.77 | 3.80 |

The initial forward round had unstable absolute timings, including tracking
20 versus 4 seconds despite unchanged graphs. A fresh unrestricted-thread serial
control then gave ID 10.98 / tracking 3.49 seconds. Therefore the variation is
**not proven to be caused by math-library threads**. Use the stable round only
as local evidence: approximately 6.5x ID speedup, but 2.7x ID+tracking batch
speedup at eight workers, excluding ingestion, saving and plots. These are laptop
measurements, not a forecast for 16/32 RCC workers or ten-year throughput.

Prefer frame parallelism with enough frames per chunk. Halo tiles provide exact
spatial seam behavior, not independently tracked subdomains. Causal tracking
remains sequential. Increasing workers cannot accelerate that state transition.

A synthetic 500,000-closed-branch stress case retained all 500,001 mappings and
the same next branch ID. Observed-only per-frame family output took about 0.0025
seconds versus 0.044 seconds enumerating the global map. This avoids repeated
historical enumeration; historical checkpoint size still grows. It is not a
measured annual-run speedup or a bounded-memory guarantee for union-find history.

## Reproduction and identity export

Run the benchmark alone, in both orders, on prepared **mm/h**, not raw RAINNC:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# macOS additionally: export VECLIB_MAXIMUM_THREADS=1
python -u scripts/benchmark_worker_modes.py PREPARED_RAIN.npy \
  --frames 24 --workers 1 2 4 8 --repeat 2 --output-dir results/workers_forward
python -u scripts/benchmark_worker_modes.py PREPARED_RAIN.npy \
  --frames 24 --workers 1 2 4 8 --repeat 2 --reverse \
  --output-dir results/workers_reverse
```

For RCC 16/32-worker tests, use at least 32 or 64 frames in a batch and compare
against a smaller-worker run on the **same node**, with the same thread limits.
Avoid simultaneous benchmarks or unrelated jobs when measuring wall time.

`branch_id` stays fixed on a continuation; split/merge creates new branches.
`family_id` is a root known at that time and can change when families merge.
Export final-known roots without changing original IDs or input catalogs:

```bash
python scripts/canonicalize_stream_families.py COMPLETED_STREAM \
  --output-dir results/final_known_families
```

The additional column is `canonical_family_id`; original `family_id` and
`branch_id` are retained. The July week export contains all 39,862 object rows.

## Freeze decision and next validation

Performance-only changes and deterministic execution have strong engineering
evidence. Keep the scientific policy opt-in: **do not promote it as a final
scientific freeze yet**. Remaining tasks are independent adjudication of
moderate-overlap continuations and ambiguous event groups, untouched historical/
future windows, and original-STEP sensitivity with consistent physical cadence
and thresholds. Avoid repeatedly tuning thresholds on these same windows.

Keep 1-mm/h identification fixed for now. Lowering to 0.1 mm/h would change
grouping and rainfall accounting, potentially increase chaining, and cannot be
justified as a tracking repair. Rain-only identification does not establish a
strong convective core. The scientific scope must be settled before a multi-year
production definition is frozen.

Tests: 102 passed, one skipped. Local NetCDF/NumPy emitted binary-compatibility
and deprecation warnings; successful reads/tests do not eliminate that environment
warning. Source hashes are stored in run settings. Tree-memory sampling lacked
process permissions locally, so reported parent RSS is not total worker memory.
