# Storm definitions: literature and operational practice

Reviewed 2026-10-02. Research only: no algorithm, preset, threshold, or tracking
graph was changed, and no new data experiment was run. This note supplements
[the original STEP anchor](ORIGINAL_STEP_ANCHOR.md) and
[our paired threshold results](RAIN_THRESHOLD_JUNE72_RESULTS.md).

## The central distinction

There is no single threshold that makes all of these objects equivalent:

| Object | Meaning in this project | What rainfall alone establishes |
|---|---|---|
| Rainfall region | A selected rain-rate support at one time | Thresholded surface precipitation |
| Detection seed / rain core | A feature used to establish identity | A relatively stronger rain feature, not certified convection |
| Rain extent | Lower-threshold pixels assigned to a seed or system | A segmentation/ownership convention |
| Branch | Unbranched temporal lineage under an explicit event convention | Algorithmic identity continuity |
| Rain event / family | Associated branches or simultaneous segments | A grouping convention, not necessarily one thunderstorm |
| Convective cell / MCS | A physically classified convective object/system | Not independently certified by RAINNC-only thresholds |

These are working distinctions, not a new universal meteorological glossary.
Our 4-km hourly RAINNC increments are interval-averaged surface precipitation.
The source contract establishes RAINNC, not any unavailable RAINC contribution.

## Verified examples, with their scope

### Original STEP: general rain events, including weak rain

Chang et al.'s author manuscript defines events through proximity, coherent
morphology and motion. It uses >0.1 mm per 3 hours (~0.033 mm/h), with 12-km
model and 4-km observations at 3-hour cadence. Large/small regions are processed
separately to limit chaining; several segments can share an event ID. Later,
events contributing the lowest collective 0.1% of rain are excluded. A strong
convective core is not mandatory. See manuscript sections 2–4:
[Chang et al., v2](https://arxiv.org/html/1601.01147).

The [publisher paper](https://journals.ametsoc.org/view/journals/clim/29/23/jcli-d-15-0844.1.xml)
was inaccessible during retrieval. Final text/supplement equivalence is not
claimed. Unit conversion does not make hourly and three-hour selection equivalent.

### Cell-TAO: deliberately selected convective rainfall objects

The reference experiment uses 5-minute data, an equivalent rate of 8.5 mm/h,
minimum area about 62 km2 and minimum lifetime 15 minutes. It varies intensity,
area, lifetime and data resolution, showing that definitions can change inferred
climate responses. This is not a threshold recommendation for hourly general
rain events. See sections 3.3–4:
[Cell-TAO, 2023](https://gmd.copernicus.org/articles/16/851/2023/).

### Operational radar: SCIT and TITAN

SCIT builds radar cells with seven reflectivity levels, 30–60 dBZ in 5-dBZ
increments, rather than treating a broad single-threshold echo as one cell.
These are cell-identification rules, not rainfall-extent thresholds. See the
[SCIT paper](https://journals.ametsoc.org/view/journals/wefo/13/2/1520-0434_1998_013_0263_tsciat_2_0_co_2.xml)
and [NWS training documentation](https://training.weather.gov/wdtd/courses/rac/documentation/rac25-products.pdf).
The NWS document was verified through its indexed excerpt; direct PDF retrieval
failed. Do not mistake historical SCIT verification criteria for all current defaults.

The inspected NCAR TITAN parameter definition defaults to 35 dBZ and minimum
size 30 km2 for 2D or 30 km3 for 3D input. Its optional dual-threshold mode finds
stronger inner regions and partitions an outer envelope. These are configurable
software settings, not a universal storm definition:
[official NCAR source](https://github.com/NCAR/lrose-core/blob/master/codebase/apps/titan/src/Titan/paramdef.Titan).
Radar reflectivity levels cannot be translated into a universal hourly rain cutoff.
The inspected repository branches are mutable, not pinned reproduction references.

### tobac: detection and extent are different tasks

tobac separates feature detection, segmentation and tracking. Multiple detection
levels resolve stronger embedded features while retaining sensitivity to weaker
ones. Watershed segmentation assigns spatial extent to seeds, potentially using
a different field. Thus a seed and its measured footprint need not have the same
threshold. This supports our architecture, not our numerical 1.0/0.5 choices:
[tobac v1.5, 2024, section 2](https://gmd.copernicus.org/articles/17/5309/2024/).

### PyFLEXTRKR: core, envelope and system classification

Cold-core detection can spread into warmer cloud extent; remaining coreless
cloud regions are also labeled. Optional coherent precipitation features combine
cloud objects. MCS classification uses tracked size/duration and precipitation
properties, retaining growth/decay stages. Bidirectional overlaps identify
merges/splits; the largest object inherits the ongoing track, with other links
recorded. Our event-ended branch IDs are a different convention:
[PyFLEXTRKR, 2023, sections 2 and 4](https://gmd.copernicus.org/articles/16/2753/2023/).

The inspected ERA5 example configuration uses PF rain threshold 0.5 mm/h and
a separate heavy-rain threshold 2 mm/h, alongside infrared cloud thresholds.
This is an application configuration, not a rain-only MCS classifier or validation
of our values:
[official ERA5 example](https://github.com/FlexTRKR/PyFLEXTRKR/blob/main/config/config_era5_mcs_tbpf.yml).

## Implications for STEP_UPGRADE — our reasoning, not source prescriptions

1. **Keep the scientific target explicit.** If studying STEP-style general
   precipitation events and their rain budget, do not silently replace that target
   with only convective cells. If selecting convection later, name and validate
   the classification separately. A 1-mm/h hourly seed should be called a rain
   seed, not a diagnosed updraft or thunderstorm core.
2. **Keep A frozen while exploring definitions.** A and B use the same tracking
   rules but different identified objects and therefore different graphs. C uses
   exactly A's graph with a lower-threshold measurement extent. A successful
   structural audit does not validate the scientific threshold.
3. **Extent expansion is defensible as a declared measurement layer.** Only
   eligible wet pixels may be assigned; dry/missing pixels are barriers. Preserve
   seeds and record ownership when multiple seeds share weak rain. Do not fill a
   circular neighborhood with invented rain or infer shared storm identity merely
   from a weak bridge. Partitioned core extent and coherent system extent answer
   different questions and should be separate outputs.
4. **The current C is stricter than a general rain-event definition.** Coreless
   weak components are unassigned, and C cannot recover weak-only birth/decay times
   outside A's graph. In the June development window it retains 87.733% of supplied
   rain versus A's 79.750% and B's 88.465%; this is coverage, not accuracy. See the
   linked local results for weighting and window limitations.
5. **Do not force a strong seed at every hour without acknowledging censoring.**
   A separate candidate could qualify a system by a substantial stronger phase
   somewhere in its lifetime while preserving associated weaker phases. This is
   a proposed experiment, not an implemented feature or a literature consensus.
6. **Hierarchy is more useful here than blindly changing the tracker.** Nested
   rain levels can expose several stronger regions inside one weak rain system.
   A hierarchy, including an HDBSCAN-style grouping proposal, still needs an
   explicit meteorological scope and temporal lineage rules. It does not resolve
   split/merge identity by itself. Investigate threshold hierarchy before adding
   an unrelated clustering dependency.
7. **Threshold and morphology interact.** Our radius-9 dilation is not a maximum
   system size or an extent distance limit. Lowering the threshold can strengthen
   chains between regions. Do not equate fewer objects with fewer false splits,
   or more persistent IDs with better tracking.

## Next bounded tests — proposed, not executed

First isolate *measurement extent* with exactly A's seeds and graph: compare
0.5 and 0.1 mm/h wet masks, preserving 1.0 as the seed/control support. Before
promoting 0.1, inspect unseeded rain, multi-seed connectivity, extreme area ratios,
rain-weighted centroid drift, dry/NaN barriers and rainfall conservation. The
existing C is uncapped in distance; large extents require review, not an arbitrary
cap presented as meteorological truth. Extent centroids must not silently drive
motion estimates.

Then, if a stronger-seed alternative is scientifically wanted, compare 1/2/5-mm/h
seed statistics independently of tracking tuning. Assess seed area and lifetime
as well as peak rate; a one-pixel crossing is not automatically a meaningful core.
Only afterwards consider changed identities or weak-stage temporal attachment.

Choose evaluation dates independently of the June development window, including
isolated cells, organized bands and broad weak rain. Compare rain coverage, event
counts, physical size, lifecycle censoring and inspected false links—not only
continuity. Use identical frozen definitions for historical/future comparisons.
Exact original-STEP comparison also needs a separately named 3-hour aggregation
control and its four-stage grouping, not just a converted threshold.

**Decision:** freeze engineering reproducibility separately from the storm
definition. Literature supports a layered seed/extent/system representation;
it does not establish 1.0, 0.5 or 0.1 mm/h as the universal correct choice.
