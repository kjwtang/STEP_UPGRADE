# STEP_UPGRADE

STEP_UPGRADE identifies two-dimensional precipitation objects and builds a
time-continuous lineage graph. Version 0.2 distinguishes three identifiers:

- a `node_id` identifies one observed object at one time;
- a `branch_id` follows an uninterrupted one-to-one track;
- a `family_id` joins branches connected by split or merge events.

The tracker supports continuation, split, merge, many-to-many complex events,
short below-threshold gaps, and checkpoint/resume across chunks and months.

## Identification

```python
import numpy as np
from step.identification import identify

precip = np.load("precip.npy")       # (time, y, x), mm per interval
valid = np.isfinite(precip)
structure = np.ones((3, 3), dtype=bool)

labels = identify(
    precip,
    morph_structure=structure,
    workers=16,
    min_size=4,
    threshold=1.0,
    valid_mask=valid,
)
```

Identification dilates a binary rain mask, labels connected regions, and
projects labels back onto observed rain pixels. The dilation footprint is a
radius around each rain pixel; it is not itself the maximum separation between
two objects. Missing pixels are barriers rather than zero precipitation.
Objects are relabeled by their first row-major rain pixel so labels do not
depend on worker scheduling.

`identify()` remains compatible with the earlier positive-rain interface. The
new `threshold` and `valid_mask` arguments make the scientific input contract
explicit. `workers` parallelizes independent frames on fork-capable HPC hosts.

## Lineage tracking

```python
from step.tracking import track_with_graph

branch_labels, graph, state = track_with_graph(
    labels,
    precip,
    tau=0.35,
    km=20,              # legacy name: grid cells per frame
    max_gap=1,          # reconnect across one missing object frame
    gap_tau=0.45,
    gap_ambiguity=0.05,
    event_overlap=0.10,
    sequence_id="present_2005_member0",
    return_state=True,
)
```

The algorithm predicts an object's next centroid from its last two observations
and uses a KD-tree to generate physically reachable candidates. Each candidate
records raw and motion-advected overlap, predicted-centroid distance, intensity
continuity, and parent/child coverage. Split and merge decisions require
overlap evidence. Remaining objects use deterministic maximum-score one-to-one
assignment.

Branch rules are intentionally strict:

| Relation | Branch behavior |
|---|---|
| one parent → one child | child keeps the branch |
| one parent → several children | parent ends; every child gets a new branch |
| several parents → one child | parents end; child gets a new branch |
| several parents → several children | parents end; every child gets a new branch |
| gap reconnection | child keeps the dormant branch |

All branches connected through an accepted split/merge event belong to the
same family. Because a later merge can join two earlier families, use the final
`graph.family_map` (or the last sequence state) as the canonical branch-to-
family mapping rather than treating an early node's provisional `family_id` as
immutable.

`track()` and the legacy arguments `phi` and `workers` remain accepted for
call compatibility. `phi` is not part of the lineage score. The legacy name
`km` still means maximum displacement in grid cells per frame; callers must
convert physical km h-1 before calling it.

## Chunk and month continuity

Never track adjacent months independently. Return state from one chunk and pass
it into the next:

```python
from step.tracking import (
    load_tracking_state,
    save_tracking_state,
    track_with_graph,
)

july_branches, july_graph, state = track_with_graph(
    july_labels, july_precip, return_state=True
)
save_tracking_state(state, "present_2005_after_july.json")

state = load_tracking_state("present_2005_after_july.json")
august_branches, august_graph, state = track_with_graph(
    august_labels,
    august_precip,
    state=state,
    return_state=True,
)
```

The checkpoint contains active and dormant object masks, motion history,
identifier counters, and family union state. It is written by atomic rename.
Tracking the same input as one block or as resumed chunks must yield identical
node IDs, branch rasters, edges, and events. Present and future climates are
separate sequences and must never share state.

For integer time coordinates, a new run starts at zero. A resumed run defaults
to `state.last_time + 1`. Use `start_time` when the first chunk must start at a
different absolute integer time.

## Validation runner

```bash
python scripts/run_npy_validation.py /path/to/precip.npy \
  --output-dir results/present_2005_07 \
  --start 0 --stop 744 \
  --threshold 1.0 --bridge-radius 9 --min-size 1 \
  --workers "${SLURM_CPUS_PER_TASK:-1}" \
  --tau 0.35 --max-displacement 20 \
  --sequence-id present_2005_member0 \
  --max-gap 1 --gap-tau 0.45 --event-overlap 0.10 \
  --gap-ambiguity 0.05 \
  --state-output results/present_2005_07/state.json

python scripts/run_npy_validation.py /path/to/precip_august.npy \
  --output-dir results/present_2005_08 \
  --threshold 1.0 --bridge-radius 9 --min-size 1 \
  --workers "${SLURM_CPUS_PER_TASK:-1}" \
  --tau 0.35 --max-displacement 20 \
  --sequence-id present_2005_member0 \
  --max-gap 1 --gap-tau 0.45 --event-overlap 0.10 \
  --gap-ambiguity 0.05 \
  --state-input results/present_2005_07/state.json \
  --state-output results/present_2005_08/state.json
```

The runner writes identified and branch rasters, node, edge and event tables,
the current family map, exact settings, and an optional next-month checkpoint.
Its input is `(time, y, x)` NumPy data. Adapt only the loading layer for
NetCDF; do not restart the tracker at a file or month boundary.

After the last month, combine the chronological tables and apply the final
canonical family map:

```bash
python scripts/finalize_catalog.py \
  results/present_2005_06 results/present_2005_07 results/present_2005_08 \
  --output-dir results/present_2005_JJA_catalog
```

The finalizer checks node uniqueness and edge references. It records the input
parts in `catalog_manifest.json`; large label rasters remain in their monthly
partitions.

## Installation and tests

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-hpc.txt
python -m pip install --no-build-isolation -e .
python -m pytest -q
```

The synthetic suite covers deterministic parallel identification, missing-data
barriers, motion, split, merge, gap expiry, JSON checkpoint round trips, and
whole-run versus resumed-run equivalence. Run the supplied Slurm smoke test on
a small real subset before production.

## Scientific validation still required

The historical upgraded identification was compared with upstream STEP on one
24-hour, 300×300 present-climate subset using a 1 mm h-1 threshold and a
9-cell dilation radius: wet-mask IoU was 1.0 and adjusted Rand index was
0.99855. This validates that subset only. Tracking, gap parameters, and event
thresholds require visual and statistical validation on real convective,
organized, weak-rain, fast-moving, present, and future cases before scientific
production.

Detailed implementation and acceptance guidance is in
[`docs/IMPLEMENTATION_GUIDE.md`](docs/IMPLEMENTATION_GUIDE.md).
