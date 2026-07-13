import numpy as np

from step.identification import identify
from step.tracking import track_with_graph


def test_parallel_identification_matches_serial():
    data = np.zeros((3, 30, 30), dtype=np.float32)
    data[:, 10:14, 10:14] = 2
    data[1, 10:14, 18:22] = 1
    structure = np.ones((3, 3), dtype=bool)
    serial = identify(data, structure, workers=1)
    parallel = identify(data, structure, workers=2)
    np.testing.assert_array_equal(serial, parallel)


def test_object_tracker_keeps_a_moving_storm_id():
    data = np.zeros((3, 40, 40), dtype=np.float32)
    for time in range(3):
        data[time, 10:15, 10 + time:15 + time] = 2
    labels = identify(data, np.ones((1, 1), dtype=bool))
    tracked, graph = track_with_graph(labels, data, tau=0.20, km=4)
    ids = [np.unique(tracked[t][tracked[t] > 0])[0] for t in range(3)]
    assert ids == [ids[0], ids[0], ids[0]]
    assert len(graph.edges) == 2
