# Storm identification/tracking: bounded literature-informed experiment

**Primary project anchor:** [Chang et al. (2016) and scaling evidence](ORIGINAL_STEP_ANCHOR.md).
The comparative methods below are secondary references. In particular, core-gated
envelopes must not silently replace STEP's general precipitation-event definition.

Date: 2026-09-26. Status: locally tested experimental implementation; no RCC
scientific validation. Scope: rain-rate identification and event-link sensitivity,
not a new convective/MCS classifier. Literature-derived observations and our design
choices are distinguished below. No external algorithm code was copied.

## 1. What other methods actually track

| Method / primary source | Object and tracking mechanism | Relevant lesson, not interchangeable thresholds |
|---|---|---|
| [TITAN, Dixon & Wiener 1993](https://journals.ametsoc.org/abstract/journals/atot/10/6/1520-0426_1993_010_0785_ttitaa_2_0_co_2.xml) | Reflectivity/size-defined contiguous radar objects; optimized correspondence, geometric split/merge treatment, motion from track history | A radar storm-cell definition is not a 1 mm/hr rain mask. History-based motion and explicit event logic are separate problems. |
| [TINT, authors' description](https://ams.confex.com/ams/98Annual/webprogram/Paper335460.html) and [official repository](https://github.com/openradar/TINT) | TITAN-inspired radar tracking; phase-correlation motion estimates and Hungarian matching | Field-based motion can supplement changing centroids. We defer implementing it until genuine hourly displacements can be tested. |
| [tobac v1.5, Sokolowsky et al. 2024](https://gmd.copernicus.org/articles/17/5309/2024/) | Modular feature detection, segmentation and temporal linking, with split/merge support | Keep seed detection, rain extent and lineage distinct; a different segmentation need not redefine every track. |
| [PyFLEXTRKR, Feng et al. 2023](https://gmd.copernicus.org/articles/16/2753/2023/) | Several applications; cloud detect-and-spread, precipitation features, overlap links in both temporal directions; explicit split/merge records | Parent/child area fractions accommodate growth and shrinkage. IoU is not the only useful overlap statistic. Rain-feature and cloud-shield thresholds are different quantities. |
| [IRT, Moseley et al. 2013](https://pure.mpg.de/rest/items/item_1850461_4/component/file_1920310/content) | Iterative rain-cell tracking estimates background flow and improves tracks iteratively | Motion estimation deserves independent validation; more complete tracks alone do not certify identity. |
| [TempestExtremes v2.1, Ullrich et al. 2021](https://gmd.copernicus.org/articles/14/5023/2021/) | General feature operators; areal tracking connects blobs through time and permits splitting/merging | Preserve explicit graph relations and do not confuse a branch with an entire connected family. |
| [Cell-TAO, Meredith et al. 2023](https://gmd.copernicus.org/articles/16/851/2023/) | Compares advection/overlap tracking and varies object area, intensity and lifetime in convection-permitting climate experiments | Their inferred climate signals can depend strongly on object-definition thresholds; preserve a fixed baseline and test sensitivity rather than tune for visual continuity. Their short-interval rates cannot be directly transplanted to hourly data. |

## 2. The 1.0 -> 0.1 mm/hr proposal

Seed-and-segment is supported as a methodological family by
[tobac 1.2](https://gmd.copernicus.org/articles/12/4551/2019/). It does **not** validate
our numerical thresholds. Here 1.0 mm/hr is an operational retained seed; 0.1 mm/hr
is eligible rain extent. Neither makes an object a deep-convective storm or MCS.
Candidate scientific questions should be separated:

1. What is the precipitation footprint associated with a retained rain seed?
2. Does this footprint belong to the same evolving rain feature at adjacent hours?
3. Does its life cycle satisfy the study's convective-storm or MCS definition?

This implementation addresses (1), and offers one controlled experiment for (2).
It does not answer (3). Stronger intensity, contiguous core area, updraft and lifetime
could inform (3), but setting values now without our data/abstract's exact criteria
would be speculative. CAPE is environmental information, not an active-storm label.

Requiring a seed every hour can truncate developing/decaying stages. Retaining
weak-only continuations needs explicit temporal hysteresis, ambiguity handling and
causal state; that is deferred. We preserve weak-only excluded rain statistics.
Do not silently call skipped weak rain “no precipitation.” No lifetime filtering is
performed on individual chunks, since that would create month-boundary artifacts.

## 3. Implemented limited changes

### E2: bounded envelope and explicit physical area

The existing partition and system experiments remain, with no default cap.
[tobac's API](https://tobac.readthedocs.io/en/stable/api/generated/tobac.segmentation.segmentation.html)
documents an optional marker-distance bound. Inspired by that control, our optional
`max_distance_km` clips eligible rain by Euclidean distance to **any retained core
pixel**, before segmentation. This is deliberately NOT a port of tobac's own-marker
distance rule, NOT a wet-path distance, and NOT a guarantee of distance from the
eventually assigned core label. Wet/dry and missing-data barriers still apply.
The cap can reduce far-away weak rain; it cannot guarantee absence of weak bridges.

Supply uniform Cartesian `(dy, dx)` spacing explicitly. Do not pass degrees or
infer 4 km from the filename. On WRF projected grids this is a nominal-distance
sensitivity experiment, not a geodesic bound. Curvilinear physical distances require
a separate implementation. Zero cap retains only seed pixels.

`min_core_area_km2` requires explicit positive scalar or 2-D cell areas. For WRF
physical area, use appropriate grid/map-factor information, not nominal dx*dy
without acknowledging the approximation. Existing `min_core_cells` remains a
separate constraint; both apply if requested. With baseline morphology, area is
summed over the grouped seed label, which may be spatially disconnected.

### E2 statistics and association fix

Save core label rasters plus per-envelope `objects.csv` and
`core_membership.csv`. Keep the full many-to-many association: a baseline seed
can contain disconnected rain patches assigned to distinct system envelopes.
It would be incorrect to assume a universal one-core-to-one-envelope relation.

Statistics include seed count, core/envelope cell counts and area ratio, rain peaks,
rain-rate sums, mean intensity, noncore rain fraction, domain-edge and missing-data
adjacency flags. These flags indicate potential censoring, not automatic rejection.
`noncore` does not necessarily mean <1 mm/hr: size-filtered intense patches can
still occur in a retained envelope. Explicit cell areas enable km2 and m3/hour;
without them physical quantities are left blank, never guessed. Volume rate is
not accumulated volume. Coreless objects are not introduced.

### T4: optional overlap-first events

Default revision-3 policy remains `score_and_overlap`. The optional `overlap`
policy removes the composite-score gate only from nontrivial event components.
It uses max(parent coverage, child coverage) and optional minimum overlap pixels.
These fractions already include raw/advected alternatives in STEP. Ordinary
one-to-one continuation and gap thresholds are unchanged.

The test runner uses overlap fraction .5 and >=2 intersecting pixels as explicit
**experimental** settings, not universal physical criteria. It first runs a matched
score-policy control with the same .5/2 settings, so comparisons isolate the score
gate. Both controls also use the same .5 candidate-overlap gate. The original .1
baseline is retained separately; changing that setting can itself change candidates.

PyFLEXTRKR motivates two-direction overlap evidence; this is not a reproduction of
that tracker. In particular, STEP ends parent branches on events rather than
inheriting the largest object's track number. Small contained fragments or advected
overlap can create spurious event edges; inspect masks and event components before
promotion. The composite score remains recorded even for accepted subthreshold
event edges. Policy is checkpointed; cross-policy resume is rejected. Default
revision-3 configuration remains checkpoint-compatible.

## 4. Deliberately deferred

* HDBSCAN replacing identification/tracking: its hierarchy can inform sensitivity
  analysis but is not itself a time-directed lineage graph.
* New optical flow, phase correlation or iterative motion estimation: would add
  another poorly constrained change before RCC validation.
* Automatically lowering tau, accepting all weak-only systems, or inventing a
  “strong core” threshold from radar or 5-minute precipitation studies.
* Automatically reclassifying whole merged families as single physical storms.
* Promoting any experimental branch to main based on synthetic tests alone.

## 5. Local evidence and limits

Synthetic tests cover contained low-score splits, small-contact rejection,
split/merge checkpoint equivalence, policy mismatch, physical area filtering,
anisotropic distance bounds, NaN/dry barriers, rain accounting, associations and
saved-data suite execution. Integration includes the pinned original STEP checkout.
An additional seeded randomized regression compared 40 eight-frame cases with
`4532608`: default rasters and complete graph records matched exactly. This is
evidence for preservation of the default path, not a real-data accuracy result.

A one-frame synthetic 1010x1634 grid with 48 seeds was exercised locally. With a
40 km cap and explicitly assumed 4 km spacing, identification plus object statistics
took about .184 s (partition) and .166 s (system); combined-process peak RSS was
about 222 MiB. This is a sparse artificial morphology, not an RCC benchmark or a
throughput forecast for a full JJA. No real-data quality claim follows from it.

## 6. Tuesday acceptance criteria

* Baseline seeds and revision-3 raster replay must exactly match saved results.
* Inspect new event edges, especially high-containment/low-IoU fragments. Count
  events, edge changes and false bridges; do not optimize for maximum count.
* Compare partition/system and .1/.2/.5 envelopes by object identity and maps,
  not color or label number. Report censored and excluded rainfall separately.
* Select any minimum area, distance and lifetime requirements in physical units,
  with justification tied to our hourly climate data and target phenomenon.
* Verify full-domain runtime and memory on RCC. No national-scale or multi-month
  performance claim until tested. Real chunk/month equivalence remains required.
* Apply identical definitions to historical and future climates; report definition
  sensitivity alongside any climate-change signal.
