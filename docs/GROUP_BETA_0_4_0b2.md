# STEP_UPGRADE group research beta: v0.4.0b2

Release date: 2026-10-05. Intended audience: internal research trials, not a public production release.
Git tag: `v0.4.0b2`; Python package version: `0.4.0b2`; distribution name: `step-upgrade`.

This is a **rainfall-object identification and lineage-tracking** beta, not a validated cloud-object
or convection classifier, and not a hosted cloud service. This release defines installation,
the recommended entry point, and delivery boundaries without changing scientific algorithms,
catalog admission criteria, or the frozen configuration. Licensing and citation preparation are
deferred for this internal trial; this is not a new licensing grant. Revisit them before external distribution.
The b2 candidate replaces the unpublished local b1 candidate with English-only repository documentation.

## 1. Recommended configuration and usage

Use the single recommended scientific preset:
`configs/rainfall_lineage_4km_1h_reference_v1.json`.

| Item | Recommended value or meaning |
|---|---|
| Grid and cadence | 4 km, hourly |
| Rainfall input | Hourly cumulative RAINNC in mm; difference snapshots to obtain mm/h |
| Catalog admission | >=1 mm/h; retain the full catalog without a stronger-core gate |
| Grouping | Morphological bridging with radius=9 grid cells; output retains actual wet pixels only |
| Tracking | Frozen overlap_ranked + coherent_adjacent configuration |
| Gap recovery | At most one frame of object absence; not a workaround for missing cumulative input files |
| 0.5/0.1 mm/h expansion | Optional separate measurement layer, not the recommended tracking feature mask |

The 9-cell radius is a grouping parameter, not a physical storm size or a maximum separation.
A label may contain disconnected rainfall fragments; it does not identify an independent convective core.
Use `run_reference_stream.py` below to load the frozen configuration explicitly.
Bare `identify()` / `track_with_graph()` defaults **do not equal** this preset;
historical experiment scripts may also use other policies.

## 2. Installation environment

Recommended Python: **3.13.5**; minimum for this beta: **3.12**.
The pinned NumPy/SciPy versions also require >=3.12. Clean-environment verification used macOS ARM64.
The execution layer uses POSIX fork and file locks; Linux/RCC users must run the smoke check on the target node.
Other Python/OS combinations have not passed a full support matrix. Native Windows is outside
the current group-beta support scope; WSL also requires independent verification.

`requirements-beta.txt` pins top-level dependencies, not every transitive dependency:
NumPy 2.5.3, SciPy 1.18.1, scikit-image 0.26.0, xarray 2026.7.0,
netCDF4 1.7.4, Matplotlib 3.11.2, and the remaining listed packages.
Save `pip freeze` for each deployment. Do not upgrade an environment during an active run or before resuming it.

### Fresh installation: macOS / Linux

Install into a new directory to avoid overwriting original STEP. `python3.13` must already be available:

```bash
git clone --branch v0.4.0b2 --single-branch \
  https://github.com/kjwtang/STEP_UPGRADE.git STEP_UPGRADE_beta_040b2
cd STEP_UPGRADE_beta_040b2
python3.13 -m venv .venv-step-beta
source .venv-step-beta/bin/activate
python -m pip install --only-binary=:all: -r requirements-beta.txt
python -m pip install --no-build-isolation -e .
python -m pip check
python -c "import step; from importlib.metadata import version; print(version('step-upgrade')); print(step.__file__)"
git describe --tags --exact-match
```

Expect package version `0.4.0b2`, tag `v0.4.0b2`, and `step.__file__` inside this checkout.
Original and upgraded STEP both use `import step`, despite different distribution names.
**Do not install original STEP and this beta into the same virtual environment.**
This editable installation requires retaining the source checkout. Do not delete it after copying
the environment, or run `git pull` / switch versions while a job is running.

If binary-wheel installation fails, record Python, CPU architecture, and the installation error.
Check wheel availability before attempting a different build; do not silently mix environments.

### RCC: existing module / conda entry point

Prepare the environment on a login node; run computation on a compute node:

```bash
module load python
source activate /project2/moyer/kjwtang/envs/wrfconda
python --version
```

Confirm Python >=3.12 first. If wrfconda is older, select a newer base Python environment
before installing these pins. Then, inside the new beta checkout:

```bash
python -m venv .venv-step-beta
source .venv-step-beta/bin/activate
python -m pip install --only-binary=:all: -r requirements-beta.txt
python -m pip install --no-build-isolation -e .
python -m pip check
```

The personal conda path above is user-specific; other group members should substitute their own base environment.
After entering a compute node, return to this checkout and reactivate `.venv-step-beta`.

## 3. First run: no research-data download required

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -u scripts/beta_smoke.py --output-dir results/beta_smoke_01 --workers 2
```

The script generates small synthetic cumulative rainfall across the January/February boundary:
12 rainfall intervals require 13 cumulative snapshots. It runs a complete control, then pauses
and resumes with a different chunk size and worker count. Complete objects, IDs, edges/scores,
events, family roots, time/grid metadata, and raster hashes must match exactly.

Expect `frames=12`, `nodes=11`, `edges=10`, `gap_edges=1`, and all eight `checks=true`.
`results/beta_smoke_01/SUCCESS` is written only after every check passes.
Use a fresh output directory for a rerun. This verifies installation and execution consistency,
not meteorological accuracy.

Full regression suite:

```bash
python -m pytest -q
```

## 4. Recommended real-data entry point

Inputs must be actual NetCDF files named, for example,
`cstm_d01_1996-07-01_00:00:00.nc`, or the same naming pattern with a `.rain` extension.
The input directory is scanned recursively. Use one climate/member sequence, with no duplicate timestamps.

Requirements:

- `RAINNC` units must be `mm`, with dimensions `(south_north, west_east)` or an additional
  singleton `Time` dimension. Only cumulative rainfall is read, not other 3D meteorological variables.
- Timestamps come from filenames. Internal `Times` and `XLAT/XLONG` arrays are not required.
- Snapshot shapes and projection attributes must agree, with `DX=DY=4000` meters. Required global attributes are
  `MAP_PROJ, CEN_LAT, CEN_LON, TRUELAT1, TRUELAT2, STAND_LON, MOAD_CEN_LAT,
  POLE_LAT, POLE_LON, DX, DY`. Do not fabricate projection metadata to pass validation.
- N hourly intervals require N+1 consecutive cumulative snapshots. Missing files, duplicate times,
  cumulative decreases/resets, and grid changes fail validation; they are not guessed, zero-filled, or repaired.
- Other rainfall components are not added automatically; RAINC, if present, is not included in RAINNC processing.

Start with a 24-hour group trial:

```bash
python -u scripts/run_reference_stream.py /PATH/TO/ONE_SEQUENCE \
  --preset configs/rainfall_lineage_4km_1h_reference_v1.json \
  --start 0 --hours 24 --sequence-id member0_1996_beta_040b2 \
  --workers 4 --chunk-frames 6 --id-executor fresh \
  --output-dir results/beta_real24_01
```

`--start` indexes the **cumulative predecessor snapshot** in chronological order,
not the first rainfall interval. The default is the full domain without cropping.
Large label arrays are not saved by default; add `--save-labels` on the initial run for small plotting cases,
and retain that option when resuming.

After interruption, repeat the command with `--resume`, keeping the same input plan, preset,
sequence ID, source, and environment. Workers, chunk size, and verified executors may change;
changing `--hours` cannot silently extend an existing plan. Do not track adjacent months independently
and splice their IDs afterward: plan a continuous interval and carry state across month boundaries.
Long plans are supported, but start trials with 24 hours; ten-year operational reliability is not established.

Progress measures processed intervals; `SUCCESS` marks completion. Multiple writers must not share an output directory.

## 5. Available capabilities

| Capability | Current status / entry point |
|---|---|
| 2D rainfall-object identification | Recommended: wet-pixel threshold, morphological grouping, deterministic labels |
| Continuous tracking, split/merge/complex events | Recommended frozen preset; produces a relationship graph, not an accuracy certification |
| Short object-absence gap recovery | Recommended preset, with ambiguity and intermediate-object conflict guards |
| Multi-worker identification | Independent frames in parallel; tracking state advances chronologically |
| Chunks, cross-month state, checkpoint/resume | Recommended runner: transactional output commits, integrity checks, single-writer lock |
| Node/branch/family IDs | Not recycled within a sequence; split/merge ends old branches and starts new ones |
| Final family catalog | `scripts/canonicalize_stream_families.py`; preserves original branch/node IDs |
| Labels, CSV, diagnostics | Runner outputs; large labels require explicit `--save-labels` |
| Exact audit across chunk/worker changes | `scripts/audit_reference_stream.py` |
| 0.5/0.1 mm/h seeded expansion | Optional measurement/research: `step.envelope.expand_seeded_envelope`; not recommended matching input |
| Stronger-rainfall phase labels | Optional offline diagnostics; no catalog/tracking change or convection classification |
| Persistent workers / tiled ID | Optional engineering experiments; check exact equivalence and performance on the target machine first |

Main chunk outputs: `objects.csv`, `edges.csv`, `events.csv`, `family_map.csv`,
`state.json`, `metadata.json`, and `COMMITTED.json` under `chunk_*/`.
Top-level outputs include `execution_contract.json`, `summary.json`, and `SUCCESS` after completion.
CSV family mappings reflect knowledge at emission time; later merges can change final family roots.
Use the final checkpoint or canonical export for complete-lifecycle analysis.

Replaying the same algorithm, inputs, and complete state should preserve IDs. IDs from runs with
different scientific configurations are not directly comparable. Filtering a display must not compact IDs.
A branch lifecycle is not an entire family lifecycle; new branch IDs after split/merge are intentional,
not necessarily tracking errors.

## 6. Unsupported or unverified claims

This beta does not validate cloud-object identification, convection classification, CAPE/updraft criteria,
cloud hosting, or automatic cumulative-reset repair. It provides no observational ground-truth accuracy,
cross-model/year generalization certification, or ten-year runtime guarantee.
Fewer breaks, more links, or more rainfall coverage do not establish higher accuracy.
Using 0.1 mm/h masks in tracking substantially changes events/gaps; it remains experimental and is not the beta default.

Start with 4 workers, 6-24-hour chunks, and short intervals. Resources depend on grid size,
object geometry, and chunk size; current evidence does not establish a universal memory requirement
or linear scaling to 16/32 workers. Checkpoint history grows with sequence length;
bounded-history cleanup is not implemented.

## 7. Verification record and feedback

See the [frozen preset](TRACKING_REFERENCE_FREEZE.md),
[continuous execution and existing JJA/winter tests](REFERENCE_STREAM_EXECUTION.md),
and [tracking sensitivity to expanded geometry](TRACKING_GEOMETRY_SENSITIVITY_2026-10-05.md).
Repeated executions are not additional independent meteorological samples.

Clean-venv installation verification on 2026-10-05 used Python 3.13.5 and macOS 15.8.1 ARM64.
Top-level dependencies were installed as binary wheels from the package index;
the modern editable build and checkout import path were verified. `pip check` reported no broken
requirements, all eight smoke exact checks passed, and the regression suite reported
**175 passed, 1 skipped** after the English-only release update, including the language regression guard.
The skipped test requires a separately pinned original STEP checkout; it is not a passed original-equivalence test.

The clean environment still emitted a NumPy/NetCDF binary-size warning and xarray/NumPy-2.5
deprecation warnings. Passing installation and tests does not eliminate native-library compatibility risk.
If a target machine shows import errors, a segmentation fault, or failed smoke checks, stop before
production calculations. Linux/RCC users must verify their own environment; macOS results are not a substitute.

Local installation and verification artifacts remain in ignored `results/beta_release_check_20261005/`,
not in the tag. The source includes the synthetic generator; no local inputs or outputs need to be copied.

Feedback should include tag/commit, `python --version`, `pip freeze`, OS/node, variable units
and cadence, preset, complete command, error logs, and a small reproducible case when available.
Do not upload raw research data, credentials, or large results to GitHub. Preserve failed outputs
for diagnosis; do not treat partial results without `SUCCESS` as a complete catalog.
