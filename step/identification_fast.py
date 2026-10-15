"""Opt-in faster equivalent of ``step.identification`` for disk footprints.

``identify_disk`` returns exactly the labels of
``identify(data, morph_structure=disk, ...)`` when the footprint is the disk
``{x^2 + y^2 <= r^2}`` used by the reference preset. Two changes, neither
altering results:

* the disk dilation is computed with an exact Euclidean distance transform,
  whose cost does not grow with the radius;
* ``binary_propagation`` is omitted: every component of the dilated, valid
  region that contains no wet cell is removed by the final ``labels * mask``,
  and propagation never changes connectivity among the remaining components.

The original module is unchanged; tests/test_identification_fast.py checks
equality on random fields with missing-data barriers.
"""

import numpy as np
from scipy import ndimage

from .identification import _stable_relabel


def disk_footprint(radius):
    yy, xx = np.ogrid[-radius:radius + 1, -radius:radius + 1]
    return xx * xx + yy * yy <= radius * radius


def identify_frame_disk(data, radius, threshold=0.0, min_size=1, valid_mask=None):
    data = np.asarray(data)
    if data.ndim != 2:
        raise ValueError("data must be two-dimensional")
    if int(radius) != radius or radius < 0:
        raise ValueError("radius must be a non-negative integer number of cells")
    if min_size < 1 or threshold < 0:
        raise ValueError("min_size must be >= 1 and threshold non-negative")
    valid = np.isfinite(data)
    if valid_mask is not None:
        valid &= np.asarray(valid_mask, dtype=bool)
    mask = valid & (data >= threshold) & (data > 0)
    if not mask.any():
        return np.zeros(mask.shape, dtype=np.int32)
    if radius == 0:
        allowed = mask
    else:
        d = ndimage.distance_transform_edt(~mask)
        allowed = (np.rint(d * d) <= radius * radius) & valid
    labels, _ = ndimage.label(allowed, structure=np.ones((3, 3), dtype=np.uint8))
    labels = _stable_relabel(labels * mask)
    if min_size > 1 and labels.max() > 0:
        counts = np.bincount(labels.ravel())
        keep = counts >= min_size
        keep[0] = False
        labels[~keep[labels]] = 0
        labels = _stable_relabel(labels)
    return labels


def identify_disk(data, radius, threshold=0.0, min_size=1, valid_mask=None):
    data = np.asarray(data)
    if data.ndim != 3:
        raise ValueError("data must have shape (time, y, x)")
    return np.stack([
        identify_frame_disk(data[t], radius, threshold, min_size,
                            None if valid_mask is None else valid_mask[t])
        for t in range(data.shape[0])
    ])
