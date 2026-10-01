#!/usr/bin/env python3
"""Bounded worker benchmark on an already prepared hourly mm/h NPY cube.

Run alone, repeat in reverse order, and constrain BLAS/OpenMP threads externally.
Does not ingest cumulative rainfall or predict RCC performance from a laptop.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
import time

import numpy as np
from step.identification import identify
from step.tracking import track_with_graph
from run_npy_validation import disk


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--frames',type=int,default=24)
    p.add_argument('--workers',type=int,nargs='+',default=[1,2,4,8])
    p.add_argument('--modes',nargs='+',choices=['frames','tiles'],default=['frames','tiles'])
    p.add_argument('--repeat',type=int,default=2)
    p.add_argument('--adjacent-policy',choices=['overlap_balanced','overlap_ranked'],default='overlap_balanced')
    p.add_argument('--reverse',action='store_true')
    a=p.parse_args()
    if not 1<=a.frames<=72 or a.repeat<1 or any(w<1 for w in a.workers):
        p.error('frames must be 1..72; repeat and workers must be positive')
    if a.output_dir.exists():raise FileExistsError(a.output_dir)
    rain=np.load(a.input,mmap_mode='r')
    if rain.ndim!=3 or len(rain)<a.frames:raise ValueError('Need a sufficient (time,y,x) cube')
    rain=rain[:a.frames]
    if np.any(np.isfinite(rain)&(rain<0)):raise ValueError('Input must already be a nonnegative hourly rate')
    a.output_dir.mkdir(parents=True)
    configs=[(m,w) for m in a.modes for w in a.workers]
    if a.reverse:configs.reverse()
    reference=None;signature=None;results=[]
    began=time.perf_counter()
    for index,(mode,workers) in enumerate(configs):
        print(f'[{index}/{len(configs)}] START {mode}, workers={workers}',flush=True)
        timings=[]
        for repeat in range(a.repeat):
            t=time.perf_counter()
            lab=identify(rain,disk(9),workers=workers,threshold=1,parallel_mode=mode)
            timings.append(time.perf_counter()-t)
            if reference is None:reference=lab.copy()
            if not np.array_equal(lab,reference):raise AssertionError('Identification differs')
            print(f'  identification repeat {repeat+1}/{a.repeat}: {timings[-1]:.3f}s',flush=True)
        t=time.perf_counter()
        tracked,graph=track_with_graph(lab,rain,adjacent_policy=a.adjacent_policy,velocity_reset='off',
            score_policy='coherent_adjacent',gap_conflict_policy='endpoint_overlap',gap_min_pixels=16,
            family_map_scope='observed')
        tracking_s=time.perf_counter()-t
        digest=hashlib.sha256(tracked.tobytes()).hexdigest()
        current=(digest,graph.objects,graph.edges,graph.events,graph.family_map)
        if signature is None:signature=current
        if current!=signature:raise AssertionError('Tracking graph or IDs differ')
        results.append(dict(mode=mode,workers=workers,identification_seconds=timings,
            median_identification_seconds=statistics.median(timings),tracking_seconds=tracking_s,
            tracked_sha256=digest,labels_graph_ids_equal=True))
        report=dict(shape=list(rain.shape),results=results,cpu_count=os.cpu_count(),adjacent_policy=a.adjacent_policy,
            thread_environment={k:os.environ.get(k) for k in ['VECLIB_MAXIMUM_THREADS',
                'OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']},
            elapsed_seconds=time.perf_counter()-began,
            note='ID and sequential tracking only; excludes ingestion/output. Local timings are not an RCC forecast.')
        (a.output_dir/'benchmark.json').write_text(json.dumps(report,indent=2))
        print(f'[{index+1}/{len(configs)}] {100*(index+1)/len(configs):.1f}% verified',flush=True)
    (a.output_dir/'SUCCESS').write_text('All worker modes reproduce reference labels, graph and IDs\n')


if __name__=='__main__':main()
