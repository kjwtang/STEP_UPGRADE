# Saved-stream adjacent pair replay

Run in the installed experimental environment:

```bash
python -u scripts/replay_stream_pairs.py \
  results/full_domain_72h_gapguard_run1 \
  --output-dir results/adjacent_replay_gapguard_run1
```

The default seven pairs concern the observed 57-61 and 60-63 cases. Pass
repeatable `--pair PARENT:CHILD` options for other node pairs. These IDs refer
to this exact saved run, not arbitrary other inputs. Requested pairs must be
adjacent in time. The saved tracking configuration is read from the final
checkpoint; there are no threshold overrides.

Reconstruct StormObject masks from saved identification rasters and use saved
weighted centroids, intensities and areas from the catalogs. This is a replay
of tracking only, NOT verification of original precipitation ingestion or
object measurement. Every output raster and each chunk's objects, edges,
events and family-map CSV must match exactly. Failure produces no diagnosis.

The report includes actual candidate admission, raw/advected overlap,
prediction error, area ratio, score, event-candidate eligibility and accepted
edges involving either endpoint. Outside-gate scores are counterfactual.
Being an event candidate alone does not establish split/merge. Rejection labels
are pair-level evidence, not a complete explanation of competing assignments.

Blocked gap candidates are annotated with `passes_gap_score`; even an
above-threshold blocked candidate is not proof it would win assignment.
Outputs are `REPORT.md` and `diagnosis.json`; sources are unchanged.
Progress counts frames replayed. Success also requires the chunk CSV checks.
