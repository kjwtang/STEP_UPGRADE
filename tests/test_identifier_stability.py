"""IDs persist within a sequence; they are not cross-run storm identities."""
import numpy as np
from step.tracking import track_with_graph, save_tracking_state, load_tracking_state


def test_expired_earlier_branch_does_not_renumber_survivor_or_reuse_id(tmp_path):
    labels = np.zeros((4, 40, 40), dtype=int)
    labels[0, 2:5, 2:5] = 1  # Earlier object disappears permanently.
    labels[0, 20:23, 20:23] = 2
    labels[1:, 20:23, 20:23] = 1  # Local labels renumber, branch must not.
    labels[3, 33:36, 33:36] = 2  # New object after the old branch expires.
    rain = (labels > 0).astype(float)*2
    first, graph, state = track_with_graph(labels[:3], rain[:3], km=2, return_state=True)
    retired, survivor = int(first[0, 3, 3]), int(first[0, 21, 21])
    assert retired != survivor
    assert np.all(first[:, 21, 21] == survivor)
    checkpoint = tmp_path/'state.json'
    save_tracking_state(state, checkpoint)
    last, _ = track_with_graph(labels[3:], rain[3:], km=2,
                               state=load_tracking_state(checkpoint))
    assert last[0, 21, 21] == survivor
    assert last[0, 34, 34] > max(retired, survivor)
    # A reporting filter selects records; it must not compact stored IDs.
    surviving_rows = [n for n in graph.objects if n.branch_id != retired]
    assert {n.branch_id for n in surviving_rows} == {survivor}
