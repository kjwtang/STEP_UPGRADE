"""Bounded, preregistered cross-date validation; does not change scientific defaults."""
from __future__ import annotations

import argparse
from datetime import datetime, timedelta, timezone
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
WINDOWS = (('july15', '1996-07-15T00:00:00'),
           ('august15', '1996-08-15T00:00:00'))


def deadline(value):
    stamp = datetime.fromisoformat(value)
    if stamp.tzinfo is None:
        raise ValueError('Deadline must include an explicit UTC offset')
    return stamp.astimezone(timezone.utc)


def commands(data, output, python, workers):
    result = []
    for name, start in WINDOWS:
        prepared, paired, extents = (output / name / child for child in
                                     ('prepared', 'paired', 'extents'))
        cube, metadata = prepared / 'rain_mm_h.npy', prepared / 'input_metadata.json'
        result.extend([
            (name, 'prepare', [python, '-u', str(ROOT / 'scripts/prepare_rain_only.py'),
                str(data), '--start', start, '--hours', '72', '--output-dir', str(prepared)]),
            (name, 'paired', [python, '-u', str(ROOT / 'scripts/compare_rain_thresholds.py'),
                str(cube), '--metadata', str(metadata), '--units', 'mm/h',
                '--workers', str(workers), '--chunk-frames', '6', '--output-dir', str(paired)]),
            (name, 'extents', [python, '-u', str(ROOT / 'scripts/compare_frozen_extents.py'),
                str(cube), '--metadata', str(metadata), '--reference', str(paired / 'A_1mm'),
                '--units', 'mm/h', '--workers', str(workers), '--chunk-frames', '6',
                '--output-dir', str(extents)]),
        ])
    return result


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def stop_group(process):
    """Terminate only this runner's newly created subprocess session."""
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()


def run(data, output, finish_before, workers=8, timeout_seconds=1200):
    if output.exists():
        raise FileExistsError(output)
    if workers < 1 or timeout_seconds <= 0:
        raise ValueError('workers and timeout must be positive')
    # Reserve time to inspect results and record the final report.
    run_until = finish_before - timedelta(minutes=10)
    if datetime.now(timezone.utc) >= run_until:
        raise ValueError('Deadline leaves no validation time plus reporting reserve')
    stages = commands(data.resolve(), output.resolve(), sys.executable, workers)
    output.mkdir(parents=True)
    dump(output / 'PLAN.json', dict(windows=[dict(name=n, predecessor=s, hours=72)
        for n, s in WINDOWS], registered_at=datetime.now(timezone.utc).isoformat(),
        finish_before=finish_before.isoformat(), subprocess_cutoff=run_until.isoformat(),
        commands=[dict(window=n, stage=s, command=c) for n, s, c in stages],
        science='Committed preset unchanged; A 1-mm seeds/graph; C 0.5/0.1 fixed-seed extents',
        scope='Held-out dates for extent threshold selection, not independent observations',
        git_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT,
                                            text=True).strip()))
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
               MKL_NUM_THREADS='1', PYTHONUNBUFFERED='1')
    statuses = []
    for index, (window, stage, command) in enumerate(stages):
        remaining = (run_until - datetime.now(timezone.utc)).total_seconds()
        if remaining <= 0:
            statuses.append(dict(window=window, stage=stage, status='deadline_not_started'))
            dump(output / 'stage_status.json', statuses)
            return False
        log = output / f'{window}_{stage}.log'
        record = dict(window=window, stage=stage, status='running', command=command,
                      log=str(log), started_at=datetime.now(timezone.utc).isoformat())
        statuses.append(record)
        dump(output / 'stage_status.json', statuses)
        print(f'START {index+1}/{len(stages)} {window}/{stage}; log={log}', flush=True)
        started = time.perf_counter()
        limit = min(timeout_seconds, remaining)
        with log.open('w') as handle:
            process = subprocess.Popen(command, cwd=ROOT, env=env, stdout=handle,
                                       stderr=subprocess.STDOUT, start_new_session=True)
            try:
                code = process.wait(timeout=limit)
                record.update(status='completed' if code == 0 else 'failed', returncode=code)
            except subprocess.TimeoutExpired:
                stop_group(process)
                record.update(status='timeout_or_deadline', returncode=process.returncode)
            except BaseException:
                stop_group(process)
                record.update(status='interrupted', returncode=process.returncode)
                record['wall_seconds'] = time.perf_counter()-started
                dump(output / 'stage_status.json', statuses)
                raise
        record['wall_seconds'] = time.perf_counter()-started
        dump(output / 'stage_status.json', statuses)
        print(json.dumps({k: v for k, v in record.items() if k != 'command'}), flush=True)
        if record['status'] != 'completed':
            return False
    # Child success contracts, not exit status alone, certify completed outputs.
    for name, _ in WINDOWS:
        for stage in ('paired', 'extents'):
            if not (output / name / stage / 'SUCCESS').is_file():
                raise AssertionError(f'Missing {name}/{stage}/SUCCESS')
    dump(output / 'SUCCESS', dict(completed_at=datetime.now(timezone.utc).isoformat(),
        scope='Execution and integrity checks passed; not scientific tracking truth'))
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data', type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--finish-before', required=True, type=deadline,
                        help='ISO timestamp including UTC offset; reserves 10 min for reporting')
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--timeout-seconds', type=float, default=1200)
    args = parser.parse_args()
    if not run(args.data, args.output_dir, args.finish_before, args.workers, args.timeout_seconds):
        raise SystemExit(2)


if __name__ == '__main__':
    main()
