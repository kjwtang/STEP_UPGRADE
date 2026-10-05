# STEP_UPGRADE group beta v0.4.0b3

Released 2026-10-05. Package: `step-upgrade==0.4.0b3`.
The frozen identification/tracking preset is unchanged from b2. This release
extends input handling, not scientific catalog admission or tracking policies.

## Installation environment

Use Python >=3.12 (reference: Python 3.13.5) on Linux/macOS in a separate environment.
The top-level dependency pins in `requirements-beta.txt` are unchanged from b2.

```bash
git clone --branch v0.4.0b3 --single-branch https://github.com/kjwtang/STEP_UPGRADE.git STEP_UPGRADE_beta_040b3
cd STEP_UPGRADE_beta_040b3
python3.13 -m venv .venv-step-beta
source .venv-step-beta/bin/activate
python -m pip install --only-binary=:all: -r requirements-beta.txt
python -m pip install --no-build-isolation -e .
python -m pip check
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -u scripts/beta_smoke.py --output-dir results/beta_smoke_01 --workers 2
```

On RCC, load the existing Python/conda environment before creating this venv:

```bash
module load python
source activate /project2/moyer/kjwtang/envs/wrfconda
```

The smoke test requires 12 frames, 11 nodes, 10 edges, one gap edge and ten true
checks. It verifies snapshot input, one multi-time file, multiple multi-time
files, a month boundary and changed worker/chunk settings on checkpoint resume.

## New input functionality

`scripts/run_reference_stream.py` accepts a directory, one NetCDF, or multiple
NetCDF paths. Default cadence is one hour; `--hours` defaults to all remaining
available intervals. Use `--inspect` to view the detected input without processing.
The README documents type detection, units, dimension overrides, file ordering,
the index-only time assumption, output schemas and restart usage.

Cumulative RAINNC is differenced across files and chunks; only the first global
snapshot lacks a predecessor. Calendar gaps, duplicate times and accumulation
resets are errors. Index-only files are assumed contiguous in file order;
this cannot establish calendar continuity independently. Supply actual grid
spacing with `--grid-km 4` if DX/DY are absent; geography is not generated.

New code/dependencies/input plans invalidate old execution contracts. Do not
update an active or paused checkout. Preserve b2 to resume b2 runs and use a
new checkout/output directory for b3.

The [b2 guide](GROUP_BETA_0_4_0b2.md) remains as a historical record.
