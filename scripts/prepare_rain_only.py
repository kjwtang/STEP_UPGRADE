#!/usr/bin/env python3
"""Prepare a bounded filename-timed RAINNC-only window for smoke testing.

Does not assert internal time or latitude/longitude validation. Original files
are read-only. Never preprocess a whole year into memory with this utility.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import xarray as xr


def prepare(directory,start,hours,output):
    if hours<1 or hours>72:
        raise ValueError('Smoke window must be 1..72 hours')
    if output.exists():
        raise FileExistsError(output)
    stamps=np.datetime64(start,'s')+np.arange(hours+1)*np.timedelta64(1,'h')
    paths=[directory/('cstm_d01_'+str(t).replace('T','_')+'.rain') for t in stamps]
    stats={p:(p.stat().st_size,p.stat().st_mtime_ns) for p in paths}
    previous=None; increments=[]; grid=None
    keys=('DX','DY','MAP_PROJ','CEN_LAT','CEN_LON','TRUELAT1','TRUELAT2','STAND_LON',
          'MOAD_CEN_LAT','POLE_LAT','POLE_LON')
    for i,path in enumerate(paths):
        with xr.open_dataset(path) as ds:
            variable=ds['RAINNC']
            if variable.dims==('Time','south_north','west_east') and variable.sizes['Time']==1:
                variable=variable.isel(Time=0)
            if variable.dims!=('south_north','west_east') or variable.attrs.get('units')!='mm':
                raise ValueError('Expected single 2-D RAINNC snapshot in mm')
            current=np.array(variable.values,dtype=np.float32,copy=True)
            signature=dict(shape=list(current.shape),**{k:float(ds.attrs[k]) for k in keys})
            if grid is not None and signature!=grid:
                raise ValueError('Grid attributes changed')
            grid=signature
        if not np.isfinite(current).all() or np.any(current<0):
            raise ValueError(f'Nonfinite or negative accumulated rain: {path.name}')
        if previous is not None:
            diff=current-previous
            if np.any(diff<0):
                raise ValueError(f'Accumulation decreased: {path.name}; min={diff.min()}, pixels={(diff<0).sum()}; no reset guessing')
            increments.append(diff)
            print(f'[{i}/{hours}] {path.name}: max increment {diff.max():.4f} mm/h',flush=True)
        previous=current
    if any((p.stat().st_size,p.stat().st_mtime_ns)!=stats[p] for p in paths):
        raise ValueError('Input changed during read; possible incomplete transfer')
    output.mkdir(parents=True)
    np.save(output/'rain_mm_h.npy',np.stack(increments))
    meta=dict(source_files=[str(p) for p in paths],source_stats=[stats[p] for p in paths],
              interval_start_times=[str(t) for t in stamps[:-1]],timestamps=[str(t) for t in stamps[1:]],
              time_source='filename only',filename_time_verified=True,cadence_hours=1,
              internal_time_verified=False,coordinates_verified=False,
              grid=grid,hours=hours,negative_increments=0,
              warning='File stability/full reads do not replace source checksums. Hourly chronology is verified against the filename contract; internal time and coordinates are not required for this grid-space test.')
    (output/'input_metadata.json').write_text(json.dumps(meta,indent=2))
    return meta


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory',type=Path)
    p.add_argument('--start',required=True,help='Timestamp of cumulative predecessor, YYYY-MM-DDTHH:MM:SS')
    p.add_argument('--hours',type=int,default=6)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args()
    prepare(a.directory,a.start,a.hours,a.output_dir)
