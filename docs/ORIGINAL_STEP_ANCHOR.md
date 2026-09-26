# Primary scientific anchor: Chang et al. (2016)

[Published paper](https://journals.ametsoc.org/view/journals/clim/29/23/jcli-d-15-0844.1.xml),
DOI 10.1175/JCLI-D-15-0844.1. Publisher full-text retrieval returned 403;
method details were checked in the [author manuscript, v2](https://arxiv.org/html/1601.01147).
The supplement and final-version method details still need verification before
claiming exact reproduction.

## Source facts (manuscript sections 2–4)

The study uses 12 km model output and 4 km Stage IV observations, with 3-hourly
analysis. The identification cutoff is >0.1 mm/3 hours (approximately .033 mm/hour),
not .1 mm/hour. Its target includes weak precipitation events, not just severe
convective cores. Identification handles large and small regions separately to
reduce chaining. Multiple simultaneous segments may share an event. The scientific
analysis separates intensity, size, duration and event count.

## Consequences for our project (design decisions, not claims from the paper)

1. This is the principal methodological baseline. Other trackers are references
   for implementation ideas, not replacements for the project's scientific target.
2. Core >=1 mm/hour plus envelope >=.1 mm/hour remains a named *alternative
   definition*. It is not a faithful recreation of the original threshold or an
   established convective classifier. Preserve a general rain-event baseline and
   account for rain excluded by every seed/extent choice.
3. Do not silently impose a strong convective-core condition on all STEP events.
   Convective classification is a separate analysis layer if required by the study.
4. Existing binary-dilation identification and sparse graph scoring are not
   numerically identical to the original four-stage identification / kernel score.
   Equal wet-mask IoU does not prove equal object grouping. Runtime improvement is
   not sufficient evidence of methodological equivalence.
5. Compare branch and family statistics separately. Neither is automatically the
   original event ID. One-to-one same-ID link counts are not directly comparable
   when the original permits several segments to share an ID.
6. Convert spatial and temporal parameters explicitly when changing grid spacing
   or sampling interval. Hourly rates and 3-hour accumulated depths are not the same
   quantity. Aggregating to the old sampling is a sensitivity control, not a reason
   to discard native resolution in production.

## Actual computational evidence

User-reported 24-frame, 600x600 RCC run, 24 GiB stage cap:

| Tracking run | Wall seconds | Parent peak RSS MiB | Interpretation |
|---|---:|---:|---|
| Original identification + original tracking | 442.13 | 10328.28 | Original tracking stage only |
| New identification + original tracking | 496.35 | 13108.21 | Same original tracker, changed grouping |
| Revision-3 tracking | 2.11 | 305.16 | Different tracker, not an exact-kernel speedup |

The first two rows exclude identification runtime; the new row is the later
revision-3 replay, not the initial .387-second core measurement. These are observed
stage measurements, not projected JJA throughput. Python reported .npy-based
tracking here, not end-to-end NetCDF loading/CAPE/3-D processing.

Inspected pinned Python reference: `ef75c083addcddbcabf6209dc95b1dff89055b0c`,
`step/tracking.py`, similarity routine. It creates a union of object coordinates,
a dense weight outer product, `squareform(pdist(union))`, and an exponential-weight
array. For U union pixels, individual float64 matrices need 8 U^2 bytes; several
coexist. At U=18147 one such matrix needs about 2.45 GiB, matching the RCC failure.
This is evidence about that Python implementation, not a claim that every
implementation of the 2016 methodology must allocate these matrices.

For a *fixed physical rain area*, a hypothetical change from 12 km to 4 km gives
roughly 9 times as many pixels, and thus 81 times as many entries in a UxU matrix.
This scaling argument assumes unchanged physical support; it is not a measured
runtime forecast. Total-domain expansion is different: a 6x domain does not imply
every object is 6x larger, nor a universal 6x or 36x runtime. Object-size distribution,
candidate counts, connected systems and I/O determine the result.

The previously inspected CSTM array is 1010x1634, about 4.58 times 600x600 cells.
The user's target domain is described as approximately 6x the test; retain that
as a planning estimate until the actual production dimensions are confirmed.

## Revised validation priorities

* **Scientific continuity:** separately report effects of grouping, similarity
  metric, event-ID convention, thresholds and time sampling. Keep historical/future
  definitions identical and check rainfall accounting.
* **Computational continuity:** test 600 crop -> full available frame -> longer
  sequence with capped memory, stage timing, max-object size and candidate counts.
  Stream 2-D rain and retain only needed tracking state. Do not load all 3-D fields.
* **Original-metric reference:** a future optimization experiment may evaluate the
  same original similarity in blocks to remove dense peak memory, before considering
  bounded-error kernel approximations. Blocking alone does not remove quadratic
  arithmetic. No such implementation or accuracy claim is added in this update.
* **No premature promotion:** T4/E2 remain optional. Their speed and synthetic
  successes do not establish equivalence to Chang et al. or full-domain validity.

This update changes the documented scientific priority, not the current algorithm.
