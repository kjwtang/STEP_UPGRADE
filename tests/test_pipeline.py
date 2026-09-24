import numpy as np
import pytest

from step.identification import identify, identify_frame
from step.tracking import load_tracking_state, save_tracking_state, track_with_graph


def labels(data):
    return identify(data, np.ones((1, 1), dtype=bool))


def test_parallel_identification_matches_serial():
    data = np.zeros((3, 30, 30), dtype=np.float32)
    data[:, 10:14, 10:14] = 2
    data[1, 10:14, 18:22] = 1
    structure = np.ones((3, 3), dtype=bool)
    np.testing.assert_array_equal(
        identify(data, structure, workers=1),
        identify(data, structure, workers=2),
    )


def test_identification_threshold_missing_barrier_and_stable_labels():
    data = np.zeros((7, 12), dtype=float)
    data[1:3, 8:10] = 2.0
    data[4:6, 1:3] = 1.0
    data[3, :] = np.nan
    result = identify_frame(data, np.ones((5, 5), dtype=bool), threshold=1.5)
    assert result.max() == 1
    assert np.all(result[1:3, 8:10] == 1)
    assert np.all(result[4:6, 1:3] == 0)
    assert np.all(result[3] == 0)


def test_object_tracker_keeps_a_moving_storm_branch():
    data = np.zeros((3, 40, 40), dtype=np.float32)
    for time in range(3):
        data[time, 10:15, 10 + time:15 + time] = 2
    tracked, graph = track_with_graph(labels(data), data, tau=0.20, km=4)
    ids = [np.unique(tracked[t][tracked[t] > 0])[0] for t in range(3)]
    assert ids == [ids[0], ids[0], ids[0]]
    assert [edge.event for edge in graph.edges] == ["continue", "continue"]


def test_split_ends_parent_and_creates_two_branches():
    data = np.zeros((2, 20, 20), dtype=float)
    data[0, 8:12, 5:15] = 2
    data[1, 8:12, 5:9] = 2
    data[1, 8:12, 11:15] = 2
    tracked, graph = track_with_graph(
        labels(data), data, tau=0.20, km=5, event_overlap=0.10
    )
    parent = np.unique(tracked[0][tracked[0] > 0])
    children = np.unique(tracked[1][tracked[1] > 0])
    assert len(parent) == 1 and len(children) == 2
    assert parent[0] not in children
    assert len(graph.events) == 1 and graph.events[0].event == "split"
    assert len([edge for edge in graph.edges if edge.event == "split"]) == 2


def test_merge_ends_parents_and_unites_families():
    data = np.zeros((2, 20, 20), dtype=float)
    data[0, 8:12, 5:9] = 2
    data[0, 8:12, 11:15] = 2
    data[1, 8:12, 5:15] = 2
    tracked, graph = track_with_graph(
        labels(data), data, tau=0.20, km=5, event_overlap=0.10
    )
    parents = np.unique(tracked[0][tracked[0] > 0])
    child = np.unique(tracked[1][tracked[1] > 0])
    assert len(parents) == 2 and len(child) == 1
    assert child[0] not in parents
    assert graph.events[0].event == "merge"
    assert len(set(graph.family_map.values())) == 1


def test_one_frame_gap_reconnects_but_longer_gap_does_not():
    data = np.zeros((6, 20, 20), dtype=float)
    data[0, 8:12, 8:12] = 2
    data[2, 8:12, 9:13] = 2
    data[5, 8:12, 10:14] = 2
    tracked, graph = track_with_graph(
        labels(data), data, tau=0.20, gap_tau=0.20, km=4, max_gap=1
    )
    first = np.unique(tracked[0][tracked[0] > 0])[0]
    second = np.unique(tracked[2][tracked[2] > 0])[0]
    third = np.unique(tracked[5][tracked[5] > 0])[0]
    assert first == second
    assert third != second
    assert [edge.event for edge in graph.edges] == ["gap_continue"]


def test_ambiguous_gap_is_rejected():
    data = np.zeros((3, 20, 20), dtype=float)
    data[0, 8:10, 4:6] = 2
    data[0, 8:10, 12:14] = 2
    data[2, 8:10, 8:10] = 2
    tracked, graph = track_with_graph(
        labels(data), data, tau=0.10, gap_tau=0.10,
        gap_ambiguity=0.05, km=10, max_gap=1,
    )
    old = set(np.unique(tracked[0][tracked[0] > 0]))
    new = np.unique(tracked[2][tracked[2] > 0])[0]
    assert new not in old
    assert not graph.edges


def test_chunk_and_checkpoint_resume_equal_whole_run(tmp_path):
    data = np.zeros((6, 24, 24), dtype=float)
    for time in (0, 1, 3, 4, 5):
        data[time, 8:12, 6 + time:10 + time] = 2
    identified = labels(data)
    whole, whole_graph = track_with_graph(
        identified, data, tau=0.20, gap_tau=0.20, km=4, max_gap=1
    )
    first, first_graph, state = track_with_graph(
        identified[:3], data[:3], tau=0.20, gap_tau=0.20, km=4,
        max_gap=1, return_state=True
    )
    checkpoint = tmp_path / "state.json"
    save_tracking_state(state, checkpoint)
    second, second_graph, _ = track_with_graph(
        identified[3:], data[3:], tau=0.20, gap_tau=0.20, km=4,
        max_gap=1, state=load_tracking_state(checkpoint), return_state=True
    )
    np.testing.assert_array_equal(whole, np.concatenate([first, second]))
    combined_edges = first_graph.edges + second_graph.edges
    assert [(e.parent_id, e.child_id, e.event) for e in whole_graph.edges] == [
        (e.parent_id, e.child_id, e.event) for e in combined_edges
    ]
    combined_nodes = first_graph.objects + second_graph.objects
    assert [(n.time, n.node_id, n.branch_id) for n in whole_graph.objects] == [
        (n.time, n.node_id, n.branch_id) for n in combined_nodes
    ]


def test_checkpoint_rejects_sequence_or_parameter_change():
    data = np.zeros((1, 10, 10), dtype=float)
    data[0, 3:6, 3:6] = 2
    identified = labels(data)
    _, _, state = track_with_graph(
        identified, data, tau=0.20, sequence_id="present_2005",
        return_state=True,
    )
    with pytest.raises(ValueError, match="sequence_id"):
        track_with_graph(
            identified, data, tau=0.20, sequence_id="future_2059",
            state=state, return_state=True,
        )
    with pytest.raises(ValueError, match="parameters"):
        track_with_graph(
            identified, data, tau=0.30, sequence_id="present_2005",
            state=state, return_state=True,
        )


def test_tracker_rejects_skipped_time_index():
    data = np.zeros((1, 10, 10), dtype=float)
    data[0, 3:6, 3:6] = 2
    identified = labels(data)
    _, _, state = track_with_graph(identified, data, return_state=True)
    with pytest.raises(ValueError, match="contiguous"):
        track_with_graph(
            identified, data, state=state, start_time=2, return_state=True
        )


def test_split_and_merge_across_chunk_boundaries_match_whole_run():
    data = np.zeros((4, 24, 24), dtype=float)
    data[0, 9:13, 5:17] = 2
    data[1, 9:13, 5:10] = 2
    data[1, 9:13, 12:17] = 2
    data[2, 9:13, 5:17] = 2
    data[3, 9:13, 6:18] = 2
    identified = labels(data)
    whole, graph = track_with_graph(
        identified, data, tau=0.20, km=5, event_overlap=0.10
    )
    chunks, events, edges, state = [], [], [], None
    for start, stop in ((0, 1), (1, 2), (2, 4)):
        raster, part, state = track_with_graph(
            identified[start:stop], data[start:stop], tau=0.20, km=5,
            event_overlap=0.10, state=state, return_state=True,
        )
        chunks.append(raster)
        events.extend(part.events)
        edges.extend(part.edges)
    np.testing.assert_array_equal(whole, np.concatenate(chunks))
    assert [(event.time, event.event) for event in graph.events] == [
        (event.time, event.event) for event in events
    ] == [(1, "split"), (2, "merge")]
    assert [(edge.parent_id, edge.child_id, edge.event) for edge in graph.edges] == [
        (edge.parent_id, edge.child_id, edge.event) for edge in edges
    ]
