"""Optional reusable frame workers; identical identification, bounded shared masks.

Create before starting other application threads. One batch completes before
the shared input buffers are refilled. No independent spatial tracking or ID
stitching is performed. Existing ``identify`` defaults remain unchanged.
"""
import multiprocessing as mp
import threading

import numpy as np

from .identification import identify, _init_pool, _identify_index

_OUTPUT = None


def _initialize_shared(rain_buffer, valid_buffer, output_buffer, shape, structure, min_size):
    global _OUTPUT
    rain = np.frombuffer(rain_buffer, dtype=np.bool_).reshape(shape)
    valid = np.frombuffer(valid_buffer, dtype=np.bool_).reshape(shape)
    rain.flags.writeable = False
    valid.flags.writeable = False
    _init_pool(rain, structure, min_size, valid)
    _OUTPUT = np.frombuffer(output_buffer, dtype=np.int32).reshape(shape)


def _identify_shared_index(index):
    _OUTPUT[index] = _identify_index(index)
    return index


class FrameIdentifier:
    """Context-managed independent-frame ID pool with reusable shared masks.

    ``capacity`` bounds the number of frames per batch; ``grid_shape`` is fixed.
    Workers receive frame indices, not serialized full rain arrays. Only wet and
    valid bool masks and disjoint output slots are shared. Workers return only
    indices, avoiding full label-array pickling. Returned labels are independent
    copies; subsequent batches cannot mutate previously returned output.
    """

    def __init__(self, grid_shape, morph_structure, capacity, workers=1,
                 threshold=0., min_size=1):
        self.grid_shape = tuple(int(v) for v in grid_shape)
        if len(self.grid_shape) != 2 or any(v < 1 for v in self.grid_shape):
            raise ValueError('grid_shape must contain two positive dimensions')
        if capacity < 1 or workers < 1 or min_size < 1 or threshold < 0:
            raise ValueError('capacity/workers/min_size must be positive; threshold nonnegative')
        self.structure = np.array(morph_structure, dtype=bool, copy=True)
        if self.structure.ndim != 2 or not self.structure.any():
            raise ValueError('morph_structure must be a non-empty 2-D array')
        self.capacity, self.workers = int(capacity), int(workers)
        self.threshold, self.min_size = threshold, min_size
        self._pool = None
        self._closed = False
        self._lock = threading.Lock()
        self._buffers = None
        if workers > 1:
            if 'fork' not in mp.get_all_start_methods():
                raise RuntimeError('Reusable frame workers require fork; use workers=1')
            context = mp.get_context('fork')
            shape = (self.capacity, *self.grid_shape)
            count = int(np.prod(shape))
            rain_buffer, valid_buffer = context.RawArray('B', count), context.RawArray('B', count)
            output_buffer = context.RawArray('i', count)
            self._buffers = (rain_buffer, valid_buffer, output_buffer)
            self._rain = np.frombuffer(rain_buffer, dtype=np.bool_).reshape(shape)
            self._valid = np.frombuffer(valid_buffer, dtype=np.bool_).reshape(shape)
            self._output = np.frombuffer(output_buffer, dtype=np.int32).reshape(shape)
            self._pool = context.Pool(processes=self.workers,
                initializer=_initialize_shared,
                initargs=(rain_buffer, valid_buffer, output_buffer, shape, self.structure, min_size))

    def identify(self, data, valid_mask=None):
        with self._lock:
            if self._closed:
                raise RuntimeError('FrameIdentifier is closed')
            data = np.asarray(data)
            if data.ndim != 3 or data.shape[1:] != self.grid_shape:
                raise ValueError('data must have the configured (time,y,x) grid')
            count = len(data)
            if not 1 <= count <= self.capacity:
                raise ValueError('Batch length must be 1..capacity')
            if valid_mask is not None and np.asarray(valid_mask).shape != data.shape:
                raise ValueError('valid_mask must have the same shape as data')
            if self._pool is None:
                return identify(data, self.structure, workers=1, threshold=self.threshold,
                                min_size=self.min_size, valid_mask=valid_mask)
            rain, valid = self._rain[:count], self._valid[:count]
            np.isfinite(data, out=valid)
            if valid_mask is not None:
                valid &= np.asarray(valid_mask, dtype=bool)
            np.greater_equal(data, self.threshold, out=rain)
            rain &= valid
            rain &= data > 0
            indices = self._pool.map(_identify_shared_index, range(count), chunksize=1)
            if indices != list(range(count)):
                raise RuntimeError('Identification workers did not complete the requested batch')
            return self._output[:count].copy()

    def close(self, terminate=False):
        with self._lock:
            if not self._closed:
                if self._pool is not None:
                    if terminate:
                        self._pool.terminate()
                    else:
                        self._pool.close()
                    self._pool.join()
                    self._pool = None
                self._buffers = None
                self._rain = self._valid = self._output = None
                self._closed = True

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.close(terminate=exc_type is not None)
