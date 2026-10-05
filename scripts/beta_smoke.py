"""Group beta installation check on synthetic cumulative rain across a month boundary."""
import argparse
from datetime import datetime, timedelta
from importlib.metadata import version
import json
from pathlib import Path
import platform
from types import SimpleNamespace

import numpy as np
import xarray as xr

from compare_rain_thresholds import reference_equivalence, dump
from run_reference_stream import DEFAULT_PRESET, run as stream_run

BETA_VERSION = '0.4.0b2'


def synthetic_input(root):
    root.mkdir()
    start = datetime(1996, 1, 31, 20)
    cumulative = np.zeros((80, 100), dtype=np.float32)
    for snapshot in range(13):
        if snapshot and snapshot != 3:  # Exactly one below-threshold object frame.
            cumulative[30:35, 20+snapshot:25+snapshot] += 2
        stamp = start+timedelta(hours=snapshot)
        attrs = {k: 1. for k in ('MAP_PROJ', 'CEN_LAT', 'CEN_LON', 'TRUELAT1',
            'TRUELAT2', 'STAND_LON', 'MOAD_CEN_LAT', 'POLE_LAT', 'POLE_LON')}
        attrs.update(DX=4000., DY=4000., note='Synthetic beta installation check, not geographic data')
        xr.Dataset({'RAINNC': (('south_north', 'west_east'), cumulative.copy(), {'units': 'mm'})},
            attrs=attrs).to_netcdf(root / f"cstm_d01_{stamp:%Y-%m-%d_%H:%M:%S}.nc")
    return root


def run(output, workers=2):
    if output.exists():
        raise FileExistsError(output)
    if workers < 1:
        raise ValueError('workers must be positive')
    installed = version('step-upgrade')
    if installed != BETA_VERSION:
        raise ValueError(f'Expected {BETA_VERSION}, installed {installed}; reinstall this beta checkout')
    output.mkdir(parents=True)
    source = synthetic_input(output / 'input')
    def options(path, **overrides):
        settings = dict(input=source, output_dir=path, preset=DEFAULT_PRESET,
            start=0, hours=12, sequence_id='synthetic_beta_month_boundary', workers=1,
            chunk_frames=12, id_executor='fresh', save_labels=True, resume=False,
            stop_after_chunks=None, no_progress=False)
        settings.update(overrides)
        return SimpleNamespace(**settings)
    print('1/3: frozen-preset whole-run control (12 hourly intervals)', flush=True)
    control = stream_run(options(output / 'control'))
    print('2/3: pause after one 4-hour chunk', flush=True)
    paused = stream_run(options(output / 'resumed', chunk_frames=4, stop_after_chunks=1))
    if paused != dict(status='paused', committed_frames=4) or (output / 'resumed/SUCCESS').exists():
        raise AssertionError('Pause was incorrectly treated as successful completion')
    print('3/3: resume with changed chunk size/workers; compare complete outputs', flush=True)
    resumed = stream_run(options(output / 'resumed', chunk_frames=3, workers=workers, resume=True))
    checks = reference_equivalence(output / 'resumed', output / 'control')
    if (control['frames'], control['nodes'], control['edges'], control['gap_edges'], control['events']) != (
            12, 11, 10, 1, {}):
        raise AssertionError(f'Synthetic reference behavior changed: {control}')
    summary = dict(status='passed', package_version=installed, python=platform.python_version(),
        platform=platform.platform(), preset=str(DEFAULT_PRESET), checks=checks,
        frames=resumed['frames'], nodes=resumed['nodes'], edges=resumed['edges'],
        gap_edges=resumed['gap_edges'], workers_on_resume=workers,
        dependency_versions={k: version(k) for k in
            ('numpy', 'scipy', 'scikit-image', 'xarray', 'netCDF4', 'pytest')},
        scope='Synthetic installation/execution check, not tracking accuracy or deployment certification')
    dump(output / 'summary.json', summary)
    dump(output / 'SUCCESS', dict(status='passed', package_version=installed))
    print(json.dumps(summary, indent=2), flush=True)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--workers', default=2, type=int)
    args = parser.parse_args()
    run(args.output_dir, args.workers)
