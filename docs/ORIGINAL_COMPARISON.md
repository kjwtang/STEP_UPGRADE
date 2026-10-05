# Original STEP vs STEP_UPGRADE: the 600x600 test domain

Historical comparison protocol. The original reference is RDCEP/STEP commit
`ef75c083addcddbcabf6209dc95b1dff89055b0c`, not an earlier STEP_UPGRADE commit.
Do not modify, install over, or replace the original STEP package. Modules are loaded under
separate names; the script verifies the reference commit and both source files against it.

## One-time preparation on a login node

Run in the separate upgraded environment:

```bash
cd /project2/moyer/kjwtang/myproject/STEP_validation/STEP_UPGRADE
git pull --ff-only
python -m pip install -r requirements-validation.txt

git clone https://github.com/RDCEP/STEP.git ../STEP_original_reference
git -C ../STEP_original_reference checkout --detach ef75c083addcddbcabf6209dc95b1dff89055b0c
```

If `STEP_original_reference` already exists, inspect its origin and status first; do not overwrite
or force-checkout it. Do not run pip install in the original directory. Original code uses the
upgraded environment's NumPy/SciPy/scikit-image. This compares original source, not a recreation
of its 2020 dependency environment. Compatibility errors remain in `run.log`; original code is not silently patched.
For the current beta, use a pinned tag/environment rather than updating a running checkout.

## Execution on a compute node

```bash
source /project2/moyer/kjwtang/myproject/STEP_validation/.venv/bin/activate
cd /project2/moyer/kjwtang/myproject/STEP_validation/STEP_UPGRADE

/usr/bin/time -v python -u scripts/compare_original_step.py \
  /project2/moyer/kjwtang/tmp/CIMP_1hr/2005.07 \
  --original-dir ../STEP_original_reference \
  --output-dir results/original_vs_upgrade_600 \
  --start 0 --hours 24 --crop 600 --threshold 1 --radius 9 \
  --workers 4 --timeout-seconds 900 --memory-gb 8
```

Scope: central 600x600 native grid, first 24 hourly increments from 25 cumulative snapshots,
1 mm/h threshold, disk radius 9. All groups share the same saved rainfall after one read/difference
pass; SHA256, source manifest, and times are recorded. Each receives identically thresholded
rainfall with min_size=1. Missing values stop comparison because original STEP lacks a missing-data
barrier; zero-filling would not be a same-input comparison. Resolution is unchanged;
the script rejects crops above 600 and full-domain mode.

| Group | Identification | Tracking | Purpose |
|---|---|---|---|
| A: original | original_id | original_track | Original end-to-end baseline |
| B: hybrid diagnosis | new_id | new_id_original_track | Same objects as C, isolating tracking differences |
| C: upgraded | new_id | new_track | Current STEP_UPGRADE |

Original tracking lacks the current checkpoint protocol, so it runs the entire 24-hour block,
not separate six-hour restarts. Upgraded tracking also processes the same whole block here.
Identification is computed once and reused, but reuse does not imply zero cost.
Core times are A=original_id+original_track, B=new_id+new_id_original_track,
and C=new_id+new_track.

## Parameters and interpretation

Original defaults `tau=0.7, phi=0.003` come from the
[original tutorial](https://github.com/relttira/STEP/wiki/Tutorial).
They are starting references, not calibrated settings for hourly 4-km data; the tutorial used
three-hour rainfall. `--original-km 30` converts 120 km / 4 km, but direction-based fallback
means it is not a strict displacement bound. Upgraded settings retain the tested
`--new-tau 0.35 --new-displacement 20` and existing gap/event settings.
Scores and motion rules differ: equal tau values would not be equivalent. This is an explicitly
configured baseline comparison, not evidence that both methods have fairly optimized parameters.
Pass calibrated `--original-tau / --original-phi / --original-km` explicitly with a new output directory
if available. Threshold and disk remain the experimental settings, not the tutorial's 0.6 or other disk shapes.

## Timeouts and memory

Each stage runs in its own subprocess with a default 900-second timeout. Five stages run serially:
up to about 75 minutes plus I/O/plotting. The 8 GiB limit combines Linux virtual-address-space
limits and sampled process-tree RSS every 0.2 seconds. Shared pages can be counted repeatedly,
so tree RSS is conservative, not exclusive physical memory. Identification workers count toward
that stage's limit. The historical suggested allocation was 4 CPUs/32 GB; ensure sufficient node time.
Outer `/usr/bin/time` MaxRSS is not the sum of exclusive memory across children.

On timeout, memory exhaustion, or error, the subprocess stops; completed results remain and
independent stages continue. Domain size and parameters do not change automatically;
missing results are not reported as zero-cost successes. Original all-pairs pixel matrices may
quickly exceed memory; this is scalability evidence, not proof of upgraded scientific accuracy.
Inspect `run.log` to distinguish MemoryError from other code/dependency errors.
`SUCCESS` means all five stages completed; `PARTIAL` indicates missing stages and exit code 2.

## Output inspection

- `identification_000.png` through `005`: shared rainfall, original grouping, upgraded grouping.
- `tracking_000.png` through `005`: rainfall and A/B/C IDs; colors/numbers are independent across columns.
- `identification_comparison.csv`: object counts, wet-mask IoU, common-wet-pixel grouping ARI; ARI is blank for empty/insufficient common pixels.
- `tracking_link_disagreements.json`: adjacent same-ID differences between B and C on identical objects.
- `tracking_summary.json`: same-ID links, shared-ID frame instances, and unique raster IDs.
- `comparison_summary.png`, `stage_results.json`: completed-stage times, failure/timeout markers, peak memory.
- Stage directories: `labels.npy`, `run.log`, parameter specifications, performance; new_track includes graph tables.
- `configuration.json`, `input_metadata.json`: versions, source/input hashes, meteorological times, provenance.

The six-frame preview corresponds to the initial plots; statistics cover all 24 frames.
Use `--plot-frames 24` for the whole period. Do not compare only ID counts or spans:
upgraded split/merge creates new branches, original tracking may share IDs across objects,
and old IDs may be reused later. Raster-ID spans are not established storm lifetimes or fragmentation rates.
Same-ID disagreement tables omit upgraded cross-branch split/merge edges; inspect
`new_track/edges.csv` as well. Check grouping in the main rainband first, then whether B/C
retain or reject connections for the same objects.
