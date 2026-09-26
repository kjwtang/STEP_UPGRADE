#!/usr/bin/env python3
"""Frame-streamed experimental envelopes on saved unthresholded rain rates."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from step.envelope import identify_envelope
from step.envelope_statistics import envelope_statistics, OBJECT_FIELDS
from run_npy_validation import disk
from validate_real_data import provenance


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('comparison',type=Path)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--core-threshold',type=float,default=1.)
    p.add_argument('--envelope-thresholds',nargs='+',type=float,default=[.1,.2,.5])
    p.add_argument('--min-core-cells',type=int,default=1)
    p.add_argument('--seed-mode',choices=['baseline','connected'],default='baseline')
    p.add_argument('--min-core-area-km2',type=float)
    areas=p.add_mutually_exclusive_group()
    areas.add_argument('--cell-area-km2',type=float,help='Explicit uniform cell area; no automatic 4 km assumption')
    areas.add_argument('--cell-area-file',type=Path,help='2-D .npy cell areas in km2 on the same crop')
    p.add_argument('--max-distance-km',type=float,help='Optional distance to ANY retained core pixel')
    p.add_argument('--spacing-km',nargs=2,type=float,metavar=('DY','DX'))
    a=p.parse_args()
    if not np.isfinite(a.core_threshold) or a.core_threshold<=0 or a.min_core_cells<1:
        p.error('Invalid core settings')
    if any(not np.isfinite(x) or not 0<x<=a.core_threshold for x in a.envelope_thresholds):
        p.error('Require 0 < envelope thresholds <= core threshold')
    source=a.comparison.resolve()
    rain=np.load(source/'rain_mm_h.npy',mmap_mode='r')
    if rain.ndim!=3 or not len(rain):
        p.error('Expected nonempty (time,y,x) rain cube')
    area=np.load(a.cell_area_file,mmap_mode='r') if a.cell_area_file else a.cell_area_km2
    config=json.loads((source/'configuration.json').read_text())
    digest=hashlib.sha256()
    for frame in rain:
        digest.update(frame.tobytes())
    if digest.hexdigest()!=config['input_sha256']:
        raise ValueError('Saved rain hash differs')
    structure=disk(config['radius']) if a.seed_mode=='baseline' else None
    options=dict(min_core_area_km2=a.min_core_area_km2,cell_area_km2=area,
                 max_distance_km=a.max_distance_km,spacing_km=a.spacing_km)
    # Validate optional physical inputs before creating output directories.
    identify_envelope(rain[0],a.core_threshold,min(a.envelope_thresholds),
                      min_core_cells=a.min_core_cells,morph_structure=structure,**options)
    a.output_dir.mkdir(parents=True,exist_ok=False)
    import scipy
    import skimage
    run_config=provenance(a)
    run_config.update(source=str(source),input_sha256=digest.hexdigest(),
                      scipy_version=scipy.__version__,skimage_version=skimage.__version__,
                      distance_definition='Euclidean distance to any retained seed pixel before segmentation')
    if a.cell_area_file:
        run_config['cell_area_sha256']=hashlib.sha256(np.asarray(area).tobytes()).hexdigest()
    (a.output_dir/'configuration.json').write_text(json.dumps(run_config,indent=2,default=str))
    saved=np.load(source/'new_id/labels.npy',mmap_mode='r')
    if saved.shape!=rain.shape:
        raise ValueError('Saved labels and rain shapes differ')
    verify=(a.seed_mode=='baseline' and a.core_threshold==config['threshold']
            and a.min_core_cells==1 and a.min_core_area_km2 is None)
    core_output=np.lib.format.open_memmap(a.output_dir/'core_labels.npy',mode='w+',dtype=np.int32,shape=rain.shape)
    # Keep labels on disk; never allocate all experiment cubes in RAM.
    started=time.perf_counter()
    for lower in sorted(set(a.envelope_thresholds)):
        for mode in ('partition','system'):
            directory=a.output_dir/f'{mode}_{lower}'; directory.mkdir()
            output=np.lib.format.open_memmap(directory/'labels.npy',mode='w+',dtype=np.int32,shape=rain.shape)
            with (directory/'statistics.csv').open('w',newline='') as f, \
                 (directory/'objects.csv').open('w',newline='') as objects_file, \
                 (directory/'core_membership.csv').open('w',newline='') as members_file:
                writer=csv.writer(f)
                object_writer=csv.DictWriter(objects_file,fieldnames=['frame']+OBJECT_FIELDS)
                member_writer=csv.DictWriter(members_file,fieldnames=['frame','envelope_label','core_label','core_cells'])
                object_writer.writeheader(); member_writer.writeheader()
                writer.writerow(['frame','core_objects','envelope_objects','core_cells',
                    'envelope_cells','unassigned_weak_cells','rain_rate_sum_assigned',
                    'rain_rate_sum_total','multi_core_envelopes'])
                for t,frame in enumerate(rain):
                    cores,env=identify_envelope(frame,a.core_threshold,lower,mode,
                        a.min_core_cells,morph_structure=structure,**options)
                    if verify and not np.array_equal(cores,saved[t]):
                        raise ValueError(f'Baseline seeds differ at frame {t}')
                    output[t]=env
                    core_output[t]=cores
                    rows,members=envelope_statistics(frame,cores,env,area)
                    object_writer.writerows(dict(frame=t,**r) for r in rows)
                    member_writer.writerows(dict(frame=t,**r) for r in members)
                    pairs=np.unique(np.column_stack((env[cores>0],cores[cores>0])),axis=0)
                    _,counts=np.unique(pairs[:,0],return_counts=True)
                    valid=np.isfinite(frame)&(frame>0)
                    writer.writerow([t,int(cores.max()),len(np.unique(env[env>0])),
                        int((cores>0).sum()),int((env>0).sum()),
                        int((valid&(frame>=lower)&(env==0)).sum()),
                        float(frame[env>0].sum(dtype=np.float64)),
                        float(frame[valid].sum(dtype=np.float64)),int((counts>1).sum())])
            output.flush(); del output
            print(f'Completed {mode} threshold={lower}',flush=True)
    core_output.flush(); del core_output
    (a.output_dir/'summary.json').write_text(json.dumps(dict(
        baseline_seeds_verified=verify,wall_seconds=time.perf_counter()-started,
        note='Experimental identification only; no tracking change. Rain-rate sums are not volume without cell areas.'),indent=2))
    (a.output_dir/'SUCCESS').touch()


if __name__=='__main__':
    main()
