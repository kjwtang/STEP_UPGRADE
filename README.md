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

## HPC validation on a `.npy` cube

Install into a fresh virtual environment, then process a small time window
before submitting a full production run:

```bash
git clone https://github.com/kjwtang/STEP_UPGRADE.git
cd STEP_UPGRADE
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-hpc.txt
python -m pip install -e .

python scripts/run_npy_validation.py /path/to/precip.npy \
  --output-dir results/smoke_24h --start 0 --stop 24 \
  --threshold 0.6 --bridge-radius 1 --min-size 4 \
  --workers "${SLURM_CPUS_PER_TASK:-1}" --tau 0.35 --max-displacement 20
```

The initial validator expects a `(time, y, x)` NumPy array. It writes
identified and tracked rasters (`.npy`), object properties (`objects.csv`),
track edges (`edges.csv`), and the exact run configuration (`run_settings.json`).

`scripts/run_validation.slurm` is a Slurm job template for the same smoke
test. Edit its resource directives and input path before submission.
