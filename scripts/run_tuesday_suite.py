#!/usr/bin/env python3
"""Bounded saved-data experiment suite; no original-tracker or NetCDF reruns."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('comparison',type=Path,help='Completed revision3 comparison directory')
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--spacing-km',type=float,nargs=2,metavar=('DY','DX'))
    p.add_argument('--max-distance-km',type=float,
                   help='Optional extra bounded-growth experiment; no default distance is assumed')
    p.add_argument('--timeout-seconds',type=float,default=900.)
    a=p.parse_args()
    import math
    if not math.isfinite(a.timeout_seconds) or a.timeout_seconds<=0:
        p.error('timeout must be finite and positive')
    if a.max_distance_km is not None and (not math.isfinite(a.max_distance_km) or
        a.max_distance_km<0 or a.spacing_km is None):
        p.error('distance cap requires a finite nonnegative value and spacing-km')
    if a.spacing_km is not None and any(not math.isfinite(x) or x<=0 for x in a.spacing_km):
        p.error('spacing must be finite and positive')
    source=a.comparison.resolve()
    if not (source/'SUCCESS').exists():
        p.error('comparison must contain SUCCESS')
    a.output_dir.mkdir(parents=True,exist_ok=False)
    scripts=Path(__file__).resolve().parent
    plans=[('events','audit_tracking_events.py',['--include-overlap-policy']),
           ('envelopes_baseline_seeds','compare_envelopes.py',[]),
           ('envelopes_connected_seeds','compare_envelopes.py',['--seed-mode','connected'])]
    if a.max_distance_km is not None:
        plans.append(('envelopes_bounded','compare_envelopes.py',
            ['--max-distance-km',str(a.max_distance_km),'--spacing-km',*map(str,a.spacing_km)]))
    statuses=[]
    for name,script,options in plans:
        command=[sys.executable,str(scripts/script),str(source),'--output-dir',
                 str(a.output_dir/name),*options]
        print('START',name,flush=True)
        started=time.perf_counter()
        with (a.output_dir/f'{name}.log').open('w') as log:
            try:
                result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,
                                      timeout=a.timeout_seconds)
                status='completed' if result.returncode==0 else 'failed'
                code=result.returncode
            except subprocess.TimeoutExpired:
                status,code='timeout',None
        statuses.append(dict(stage=name,status=status,returncode=code,command=command,
                             wall_seconds=time.perf_counter()-started))
        (a.output_dir/'suite_status.json').write_text(json.dumps(statuses,indent=2))
        if status!='completed':
            raise SystemExit(f'{name}: {status}; inspect its log, no suite SUCCESS written')
    (a.output_dir/'SUCCESS').touch()
    print('Completed:',a.output_dir,flush=True)


if __name__=='__main__':
    main()
