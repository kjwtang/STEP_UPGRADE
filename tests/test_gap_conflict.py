import numpy as np
import pytest
from step.tracking import track_with_graph, save_tracking_state, load_tracking_state


def example(middle=True):
    lab = np.zeros((3, 20, 20), dtype=int)
    lab[0, 5:10, 5:10] = lab[2, 5:10, 5:10] = 1
    if middle:
        lab[1, 2:15, 2:15] = 1
    return lab, (lab > 0).astype(float)*2


def test_guard_blocks_observed_middle_but_preserves_true_gap():
    for middle in (True, False):
        lab, rain = example(middle)
        _, base = track_with_graph(lab, rain, tau=.99, gap_tau=.45)
        assert any(e.gap for e in base.edges)
        conflicts = []
        _, graph = track_with_graph(lab, rain, tau=.99, gap_tau=.45,
            gap_conflict_policy='endpoint_overlap', gap_conflict_callback=conflicts.extend)
        assert bool(conflicts) == middle
        assert any(e.gap for e in graph.edges) == (not middle)


def test_guard_checkpoint_equivalence_and_policy_mismatch(tmp_path):
    lab, rain = example()
    opts = dict(tau=.99, gap_tau=.45, gap_conflict_policy='endpoint_overlap')
    expected, graph = track_with_graph(lab, rain, **opts)
    first, g1, state = track_with_graph(lab[:2], rain[:2], return_state=True, **opts)
    path = tmp_path/'state.json'
    save_tracking_state(state, path)
    restored = load_tracking_state(path)
    last, g2 = track_with_graph(lab[2:], rain[2:], state=restored, **opts)
    np.testing.assert_array_equal(expected, np.concatenate([first, last]))
    assert graph.edges == g1.edges + g2.edges
    with pytest.raises(ValueError, match='parameters differ'):
        track_with_graph(lab[2:], rain[2:], state=load_tracking_state(path), tau=.99)


def test_restrict_guard_to_one_missing_frame():
    lab, rain = example()
    with pytest.raises(ValueError, match='max_gap=1'):
        track_with_graph(lab, rain, max_gap=2, gap_conflict_policy='endpoint_overlap')
