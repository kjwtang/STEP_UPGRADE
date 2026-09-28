# Streaming frame progress

`validate_real_data.py --stream` now shows an overall frame bar after each tracked
frame, accumulated across chunk boundaries. Hourly input uses the label `hours`;
other cadence uses `frames`. Example (illustrative elapsed time):

```
[############------------] 12/24 hours tracked (50.0%) | elapsed 35s | ETA ~35s | tracking
```

Reading, identification, checkpoint saving and report generation are named phases.
The percentage is completed **tracking frames**, not a weighted estimate of total
CPU work or durable output. At 100%, final saving/reporting may still be running.
Only `SUCCESS: outputs saved` confirms the SUCCESS marker was written. On failure,
the counter remains at its last value and is never forced to 100%.

ETA is an approximate extrapolation from elapsed time and completed frames; it
includes preceding I/O but can fluctuate with storm density, chunks and plots.
Progress is printed to stderr as log-friendly lines. No extra package required.
`--no-progress` disables it. This update affects the streaming validation command,
not the separate multi-experiment suite runner.

RCC result reported by user: 24 frames, 1010x1634, default tracking, six-frame
chunks, four identification workers. External wall time 74.29 s; reported stage
totals: read 26.55 s, identification 4.37 s, tracking 29.71 s, output 2.71 s,
preview plots 9.19 s. Sampled summed tree RSS 1680.23 MiB includes shared fork
pages more than once. This run had 3266 object nodes, 738 edges including 3 gap
edges, and no split/merge events. Whole-versus-chunk equivalence was not checked
in this streaming run. These are results of one day, not a full-season forecast.
