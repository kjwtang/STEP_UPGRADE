# One-frame gap conflict experiment

This is an opt-in diagnostic intervention, not a validated scientific default.
Identification, adjacent matching, tau, and event policy are unchanged.

The 72-hour RCC audit showed endpoint gaps across intermediate identified
objects, including alternating 60->62 and 61->63 chains. The saved audit only
measures overlap with the union of endpoint footprints. It does not establish
why individual adjacent edges failed or whether split/merge was appropriate.

`--gap-conflict-policy endpoint_overlap` blocks a candidate gap when the SAME
immediately preceding object intersects at least 50% AND 16 pixels of EACH
endpoint separately, without advection. These are experimental guard constants,
not literature-calibrated storm thresholds. The default is `off`. Only
max_gap=1 is supported, to avoid pretending that older intermediate frames are
available in the checkpoint. Small endpoints cannot trigger this guard.

Stream validation writes `chunk_*/gap_conflicts.json`: these are blocked
CANDIDATES, not necessarily edges baseline assignment would have selected.
The guard does not create replacement edges. Fewer gaps can mean more broken
tracks, not better tracking. The policy is stored in checkpoint configuration;
cross-policy resumption is rejected. Existing baseline checkpoints remain valid
for baseline runs.

Compare the same 72-hour input with baseline settings, changing only the policy
and output directory. Check runtime, gap/edge counts, conflict records and the
57-61, 60-63, 65-68 case maps. Do not promote this policy merely because gap count
falls. Next investigate rejected adjacent candidate scores/gates and morphology
changes, then evaluate replacements with split/merge-aware tests.

Local tests cover true empty-frame gaps, observed intermediate conflicts,
checkpoint raster/edge equivalence, and cross-policy checkpoint rejection.
RCC scientific validation is pending. Images use identified masks, not rain
intensity; no inference of actual missing observation is warranted.
