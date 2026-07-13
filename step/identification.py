"""Parallel, two-dimensional precipitation-object identification.

The public ``identify`` signature remains compatible with STEP.  Unlike the
original implementation, morphology is applied to a binary rain mask rather
than to an integer label image, so results cannot depend on label values.
"""

import multiprocessing as mp

import numpy as np
from scipy import ndimage
from skimage.segmentation import relabel_sequential


_DATA = None
_STRUCTURE = None
_MIN_SIZE = None


def _init_pool(data, structure, min_size):
    global _DATA, _STRUCTURE, _MIN_SIZE
    _DATA = data
    _STRUCTURE = structure
    _MIN_SIZE = min_size


def _identify_slice(mask, structure, min_size):
    """Join nearby rain pixels while retaining only pixels observed as rain."""
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return np.zeros(mask.shape, dtype=np.int32)

    # Dilating a *binary* mask defines a physically configurable joining gap.
    # Reapplying the original mask prevents artificial rain pixels in output.
    joined = ndimage.binary_dilation(mask, structure=structure)
    labels, _ = ndimage.label(joined, structure=np.ones((3, 3), dtype=np.uint8))
    labels = labels * mask
    labels = relabel_sequential(labels)[0].astype(np.int32, copy=False)

    if min_size > 1 and labels.max() > 0:
        counts = np.bincount(labels.ravel())
        keep = counts >= min_size
        keep[0] = False
        labels[~keep[labels]] = 0
        labels = relabel_sequential(labels)[0].astype(np.int32, copy=False)
    return labels


def _identify_index(index):
    return _identify_slice(_DATA[index] > 0, _STRUCTURE, _MIN_SIZE)


def identify(data, morph_structure=None, workers=1, min_size=1, chunksize=1):
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

    structure = np.asarray(morph_structure if morph_structure is not None else np.ones((3, 3)), dtype=bool)
    if structure.ndim != 2 or not structure.any():
        raise ValueError("morph_structure must be a non-empty 2-D array")

    count = data.shape[0]
    if workers == 1 or count < 2:
        frames = [_identify_slice(data[i] > 0, structure, min_size) for i in range(count)]
    else:
        methods = mp.get_all_start_methods()
        if "fork" not in methods:
            raise RuntimeError("parallel identify requires a fork-capable platform; use workers=1")
        context = mp.get_context("fork")
        with context.Pool(processes=workers, initializer=_init_pool, initargs=(data, structure, min_size)) as pool:
            frames = pool.map(_identify_index, range(count), chunksize=chunksize)
    return np.stack(frames, axis=0)
