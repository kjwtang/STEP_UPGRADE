#!/usr/bin/env python3
"""Controlled fresh/persistent ID batch benchmark; no ingestion or tracking speed claim."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time

import numpy as np
from step.identification import identify
from step.parallel_identification import FrameIdentifier
from run_npy_validation import disk


def hashes(labels):
    return [hashlib.sha256(memoryview(frame)).hexdigest() for frame in labels]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input', type=Path, help='Already normalized hourly mm/h NPY cube')
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--frames', type=int, default=24)
    p.add_argument('--batch-frames', type=int, nargs='+', default=[6, 24])
    p.add_argument('--workers', type=int, nargs='+', default=[4, 8])
    p.add_argument('--repeat', type=int, default=2)
    p.add_argument('--reverse', action='store_true')
    a = p.parse_args()
    if not 1 <= a.frames <= 72 or a.repeat < 1 or min(a.workers + a.batch_frames) < 1:
        p.error('frames must be 1..72; workers/batch/repeat positive')
    source = np.load(a.input, mmap_mode='r')
    if source.ndim != 3 or len(source) < a.frames:
        raise ValueError('Need sufficient prepared (time,y,x) frames')
    rain = source[:a.frames]
    if np.any(np.isfinite(rain) & (rain < 0)):
        raise ValueError('Already normalized nonnegative rain is required')
    a.output_dir.mkdir(parents=True, exist_ok=False)
    expected = []
    # Serial preflight supplies exact label hashes and warms the same prepared
    # pages; retain only hashes, not a second complete cube of labels.
    for offset in range(0, a.frames, 6):
        expected.extend(hashes(identify(rain[offset:offset + 6], disk(9), threshold=1)))
    configs = [(batch, workers, executor) for batch in a.batch_frames
               for workers in a.workers for executor in ('fresh', 'persistent')]
    if a.reverse:
        configs.reverse()
    records = []
    for index, (batch, workers, executor) in enumerate(configs):
        for repeat in range(a.repeat):
            print(f'[{index + 1}/{len(configs)}] batch={batch}, workers={workers}, '
                  f'{executor}, repeat={repeat + 1}/{a.repeat}', flush=True)
            service = None
            t = time.perf_counter()
            if executor == 'persistent':
                service = FrameIdentifier(rain.shape[1:], disk(9), capacity=batch,
                                          workers=workers, threshold=1)
            startup = time.perf_counter() - t
            per_batch, actual = [], []
            try:
                for offset in range(0, a.frames, batch):
                    data = rain[offset:offset + batch]
                    t = time.perf_counter()
                    labels = (identify(data, disk(9), workers=workers, threshold=1)
                              if service is None else service.identify(data))
                    per_batch.append(time.perf_counter() - t)
                    actual.extend(hashes(labels))
                    del labels
                if actual != expected:
                    raise AssertionError('Frame labels differ from serial reference')
            finally:
                t = time.perf_counter()
                if service is not None:
                    service.close()
                shutdown = time.perf_counter() - t
            records.append(dict(batch_frames=batch, workers=workers, executor=executor,
                repeat=repeat, startup_seconds=startup, shutdown_seconds=shutdown,
                batch_identification_seconds=per_batch,
                total_identification_seconds=startup + sum(per_batch) + shutdown,
                per_frame_sha256_equal=True))
            report = dict(shape=list(rain.shape), records=records, reverse=a.reverse,
                cpu_count=os.cpu_count(),
                thread_environment={key: os.environ.get(key) for key in
                    ('VECLIB_MAXIMUM_THREADS', 'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS')},
                note='Prepared warm-cache ID only. Startup and shutdown included; digest, ingestion, tracking and output excluded. Local laptop timing is not RCC scaling evidence.')
            (a.output_dir / 'benchmark.json').write_text(json.dumps(report, indent=2))
    (a.output_dir / 'SUCCESS').write_text('All frame hashes equal serial reference\n')


if __name__ == '__main__':
    main()
