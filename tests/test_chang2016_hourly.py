import numpy as np
from scipy import ndimage

from step import chang2016 as c
from step.chang2016_hourly import mask_distance_cells, track_hourly


def _blob(shape, cy, cx, r, value=2.0):
    yy, xx = np.ogrid[:shape[0], :shape[1]]
    return np.where((yy - cy) ** 2 + (xx - cx) ** 2 <= r * r, value, 0.0)


def _labels(rain):
    out = np.zeros(rain.shape, dtype=np.int32)
    for t in range(len(rain)):
        out[t] = ndimage.label(rain[t] > 0, structure=np.ones((3, 3)))[0]
    return out


def _run(frames, **kwargs):
    rain = np.stack(frames)
    return track_hourly(_labels(rain), rain, grid_km=4, **kwargs)


def _segments(rain):
    return c.segments_from_labels(ndimage.label(rain > 0, structure=np.ones((3, 3)))[0], rain)


def test_mask_distance_matches_bruteforce():
    rng = np.random.default_rng(3)
    shape = (60, 90)
    for _ in range(5):
        a = _blob(shape, rng.integers(10, 50), rng.integers(10, 40), rng.integers(2, 8))
        b = _blob(shape, rng.integers(10, 50), rng.integers(50, 80), rng.integers(2, 8))
        sa, sb = _segments(a)[0], _segments(b)[0]
        dy = sa.ys[:, None] - sb.ys[None, :]
        dx = sa.xs[:, None] - sb.xs[None, :]
        expected = float(np.sqrt(dy * dy + dx * dx).min())
        assert mask_distance_cells(sa, sb, 1000.0) == expected
        got = mask_distance_cells(sa, sb, expected - 0.5)
        assert np.isinf(got)


def test_moving_storm_keeps_one_event():
    shape = (80, 160)
    events, rows, _ = _run([_blob(shape, 40, 20 + 3 * t, 6) for t in range(5)])
    assert len({r["event_id"] for r in rows}) == 1
    assert all(r["mask_distance_km"] == 0 for r in rows if r["rule"] != "new")


def test_far_segment_is_a_new_event():
    shape = (80, 200)
    frames = [_blob(shape, 40, 20, 6), _blob(shape, 40, 20 + 40, 6)]  # edges 28 cells (112 km) apart
    _, rows, _ = _run(frames)
    assert len({r["event_id"] for r in rows}) == 2


def test_gap_link_keeps_event_and_is_optional():
    shape = (60, 80)
    frames = [_blob(shape, 30, 30, 6), np.zeros(shape), _blob(shape, 30, 32, 6)]
    _, rows, _ = _run(frames)
    assert len({r["event_id"] for r in rows}) == 1
    assert [r["gap"] for r in rows] == [0, 1]
    _, rows, _ = _run(frames, allow_gap=False)
    assert len({r["event_id"] for r in rows}) == 2


def test_split_shares_event():
    shape = (80, 120)
    t0 = _blob(shape, 40, 60, 10)
    t1 = _blob(shape, 25, 60, 6) + _blob(shape, 55, 60, 6)
    _, rows, _ = _run([t0, t1], max_distance_km=60.0)
    assert len({r["event_id"] for r in rows}) == 1
    assert sum(r["time"] == 1 for r in rows) == 2


def test_most_similar_parent_wins_over_largest():
    shape = (80, 140)
    big, small = _blob(shape, 40, 30, 12), _blob(shape, 40, 80, 5)
    child = _blob(shape, 40, 70, 5)
    events, rows, _ = _run([big + small, child], max_distance_km=140.0)
    assert events[1, 40, 70] == events[0, 40, 80] != events[0, 40, 30]


def test_rows_feed_factorize():
    shape = (80, 160)
    _, rows, _ = _run([_blob(shape, 40, 20 + 3 * t, 6) for t in range(4)])
    f = c.factorize(rows, 16.0, 1.0, 4)
    assert f["number_of_events"] == 1 and f["duration_factor_h"] == 4.0
