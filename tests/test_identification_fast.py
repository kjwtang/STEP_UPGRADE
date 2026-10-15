import numpy as np
import pytest
from scipy import ndimage

from step.identification import identify
from step.identification_fast import disk_footprint, identify_disk


@pytest.mark.parametrize("radius", [0, 1, 3, 9])
@pytest.mark.parametrize("min_size", [1, 5])
def test_identify_disk_equals_reference_identify(radius, min_size):
    rng = np.random.default_rng(radius + 10 * min_size)
    frames = []
    for _ in range(4):
        f = ndimage.gaussian_filter(rng.standard_normal((150, 170)), 4) * 10
        frames.append(f)
    data = np.stack(frames)
    data[0, :, 80] = np.nan          # missing-data barrier
    valid = rng.random(data.shape) > 0.01
    valid[1, 40:44, :] = False
    want = identify(data, morph_structure=disk_footprint(radius), threshold=1.0,
                    min_size=min_size, valid_mask=valid)
    got = identify_disk(data, radius, threshold=1.0, min_size=min_size, valid_mask=valid)
    assert np.array_equal(want, got)
