import numpy as np
import pytest

from step.identification import identify
from step.parallel_identification import FrameIdentifier


@pytest.mark.parametrize('workers', [1, 2, 4])
@pytest.mark.parametrize('min_size', [1, 4])
def test_repeated_shared_batches_match_serial(workers, min_size):
    rng = np.random.default_rng(992)
    footprint = np.ones((4, 5), dtype=bool)  # even footprint is intentional
    with FrameIdentifier((23, 31), footprint, capacity=4, workers=workers,
                         threshold=1., min_size=min_size) as service:
        previous = None
        for count in (4, 1, 3, 2, 4):
            rain = rng.uniform(0, 2, (count, 23, 31))
            rain[rng.random(rain.shape) < .1] = np.nan
            valid = rng.random(rain.shape) > .05
            expected = identify(rain, footprint, threshold=1., min_size=min_size, valid_mask=valid)
            actual = service.identify(rain, valid_mask=valid)
            np.testing.assert_array_equal(actual, expected)
            if previous is not None:
                np.testing.assert_array_equal(previous[0], previous[1])
            previous = (actual, actual.copy())
    with pytest.raises(RuntimeError, match='closed'):
        service.identify(rain)


def test_no_stale_wet_masks_after_empty_and_missing_frames():
    with FrameIdentifier((12, 13), np.ones((3, 3)), 3, workers=2, threshold=1) as service:
        assert service.identify(np.full((3, 12, 13), 2.)).max() == 1
        assert not service.identify(np.zeros((1, 12, 13))).any()
        assert not service.identify(np.full((2, 12, 13), np.nan)).any()
        assert not service.identify(np.full((3, 12, 13), .5)).any()
        assert service.identify(np.full((1, 12, 13), 2.)).max() == 1


def test_worker_pool_is_reused_and_capacity_is_enforced():
    with FrameIdentifier((8, 9), np.ones((3, 3)), 2, workers=2) as service:
        pool = service._pool
        pids = [worker.pid for worker in pool._pool]
        service.identify(np.ones((1, 8, 9)))
        service.identify(np.ones((2, 8, 9)))
        assert [worker.pid for worker in pool._pool] == pids
        with pytest.raises(ValueError, match='length'):
            service.identify(np.ones((3, 8, 9)))
        with pytest.raises(ValueError, match='grid'):
            service.identify(np.ones((1, 9, 8)))
        with pytest.raises(ValueError, match='valid_mask'):
            service.identify(np.ones((1, 8, 9)), valid_mask=np.ones((8, 9)))
    service.close()  # idempotent
