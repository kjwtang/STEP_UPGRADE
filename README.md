# STEP_UPGRADE

Scalable replacement for the slowest parts of STEP: binary 2-D segmentation
and sparse object-based temporal tracking.  It preserves the principal STEP
function signatures while avoiding pixel-by-pixel all-pairs distances.

## HPC use

```python
from step.identification import identify
from step.tracking import track_with_graph

labels = identify(precip, morph_structure, workers=16, min_size=4)
tracked, graph = track_with_graph(labels, precip, tau=0.35, km=20)
```

`workers` parallelizes timestep segmentation with Linux `fork`, which is the
expected HPC environment.  Set it no higher than the CPUs allocated by Slurm
(`SLURM_CPUS_PER_TASK`). Tracking is time-ordered and intentionally remains
serial; its candidate comparisons are sparse and avoid STEP's all-pixel
distance calculation.

## Important compatibility notes

- Input is `(time, y, x)` and positive precipitation values are rain.
- `morph_structure` is a 2-D binary gap-bridging structure, as in STEP.
- `track(..., phi=...)` accepts `phi` for legacy callers but does not use it.
- The accompanying `TrackGraph` preserves merge edges. Future versions will
  add split and short-gap reactivation edges.

## Test

```bash
python -m pytest -q
```
