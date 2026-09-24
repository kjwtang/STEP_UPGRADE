"""Parallel, two-dimensional precipitation-object identification.

The public ``identify`` signature remains compatible with STEP.  Unlike the
original implementation, morphology is applied to a binary rain mask rather
than to an integer label image, so results cannot depend on label values.
"""

import multiprocessing as mp

import numpy as np
from scipy import ndimage


_DATA = None
_STRUCTURE = None
_MIN_SIZE = None
_VALID = None


def _init_pool(data, structure, min_size, valid):
    global _DATA, _STRUCTURE, _MIN_SIZE, _VALID
    _DATA = data
    _STRUCTURE = structure
    _MIN_SIZE = min_size
    _VALID = valid


def _stable_relabel(labels):
    """Relabel objects by their first row-major rain pixel."""
    present = np.unique(labels)
    present = present[present > 0]
    if not len(present):
        return np.zeros(labels.shape, dtype=np.int32)
    ordered = sorted(
        present.tolist(),
        key=lambda value: int(np.flatnonzero(labels == value)[0]),
    )
    lookup = np.zeros(int(labels.max()) + 1, dtype=np.int32)
    for new, old in enumerate(ordered, start=1):
        lookup[old] = new
    return lookup[labels]


def _identify_slice(mask, structure, min_size, valid_mask=None):
    """Join nearby rain pixels while retaining only pixels observed as rain."""
    mask = np.asarray(mask, dtype=bool)
    valid = (
        np.ones(mask.shape, dtype=bool)
        if valid_mask is None
        else np.asarray(valid_mask, dtype=bool)
    )
    if valid.shape != mask.shape:
        raise ValueError("valid_mask must have the same shape as the data")
    mask &= valid
    if not mask.any():
        return np.zeros(mask.shape, dtype=np.int32)

    # Dilating a *binary* mask defines a physically configurable joining gap.
    # Reapplying the original mask prevents artificial rain pixels in output.
    allowed = ndimage.binary_dilation(mask, structure=structure) & valid
    # Propagation through the valid part of the dilated region prevents a
    # footprint from jumping across a missing-data barrier.
    joined = ndimage.binary_propagation(
        mask,
        structure=np.ones((3, 3), dtype=bool),
        mask=allowed,
    )
    labels, _ = ndimage.label(joined, structure=np.ones((3, 3), dtype=np.uint8))
    labels = labels * mask
    labels = _stable_relabel(labels)

    if min_size > 1 and labels.max() > 0:
        counts = np.bincount(labels.ravel())
        keep = counts >= min_size
        keep[0] = False
        labels[~keep[labels]] = 0
        labels = _stable_relabel(labels)
    return labels


def _identify_index(index):
    valid = None if _VALID is None else _VALID[index]
    return _identify_slice(_DATA[index] > 0, _STRUCTURE, _MIN_SIZE, valid)


def identify(
    data,
    morph_structure=None,
    workers=1,
    min_size=1,
    chunksize=1,
    threshold=0.0,
    valid_mask=None,
):
    """Identify 2-D precipitation objects independently at each timestep.

    Parameters are deliberately in grid cells for backward compatibility.
    ``morph_structure`` is normally a disk or square representing the maximum
    gap to bridge.  On Linux/HPC, ``workers > 1`` uses forked worker processes;
    the input cube is inherited copy-on-write and is not serialized per frame.
    """
    data = np.asarray(data)
    if data.ndim != 3:
        raise ValueError("data must have shape (time, y, x)")
    if workers < 1:
        raise ValueError("workers must be at least one")
    if min_size < 1:
        raise ValueError("min_size must be at least one")
    if threshold < 0:
        raise ValueError("threshold must be non-negative")

    valid = np.isfinite(data)
    if valid_mask is not None:
        supplied = np.asarray(valid_mask, dtype=bool)
        if supplied.shape != data.shape:
            raise ValueError("valid_mask must have the same shape as data")
        valid &= supplied
    rain = valid & (data >= threshold) & (data > 0)

    structure = np.asarray(morph_structure if morph_structure is not None else np.ones((3, 3)), dtype=bool)
    if structure.ndim != 2 or not structure.any():
        raise ValueError("morph_structure must be a non-empty 2-D array")

    count = data.shape[0]
    if workers == 1 or count < 2:
        frames = [
            _identify_slice(rain[i], structure, min_size, valid[i])
            for i in range(count)
        ]
    else:
        methods = mp.get_all_start_methods()
        if "fork" not in methods:
            raise RuntimeError("parallel identify requires a fork-capable platform; use workers=1")
        context = mp.get_context("fork")
        with context.Pool(
            processes=workers,
            initializer=_init_pool,
            initargs=(rain, structure, min_size, valid),
        ) as pool:
            frames = pool.map(_identify_index, range(count), chunksize=chunksize)
    return np.stack(frames, axis=0)


def identify_frame(data, morph_structure=None, min_size=1, threshold=0.0,
                   valid_mask=None):
    """Identify one 2-D frame with explicit threshold and validity."""
    data = np.asarray(data)
    if data.ndim != 2:
        raise ValueError("data must be two-dimensional")
    return identify(
        data[np.newaxis],
        morph_structure=morph_structure,
        workers=1,
        min_size=min_size,
        threshold=threshold,
        valid_mask=(
            None if valid_mask is None else np.asarray(valid_mask)[np.newaxis]
        ),
    )[0]
