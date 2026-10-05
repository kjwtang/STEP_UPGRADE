from datetime import datetime, timedelta, timezone
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import validate_extent_holdouts as runner


def test_preregistered_windows_and_unchanged_commands(tmp_path):
    stages = runner.commands(tmp_path / 'data', tmp_path / 'out', 'python', 8)
    assert len(stages) == 6
    assert [c[c.index('--start')+1] for _, s, c in stages if s == 'prepare'] == [
        '1996-07-15T00:00:00', '1996-08-15T00:00:00']
    for _, stage, command in stages:
        if stage != 'prepare':
            assert command[command.index('--workers')+1] == '8'
            assert command[command.index('--chunk-frames')+1] == '6'
            assert '--preset' not in command and '--threshold' not in command
    assert runner.deadline('2026-10-05T14:00:00-05:00').hour == 19
    with pytest.raises(ValueError, match='offset'):
        runner.deadline('2026-10-05T14:00:00')


def test_guarded_deadline_and_existing_output(tmp_path):
    past = datetime.now(timezone.utc)
    with pytest.raises(ValueError, match='Deadline'):
        runner.run(tmp_path, tmp_path / 'none', past)
    assert not (tmp_path / 'none').exists()
    with pytest.raises(FileExistsError):
        runner.run(tmp_path, tmp_path, past)


@pytest.mark.parametrize('code,timeout,status', [
    ('raise SystemExit(4)', 5, 'failed'),
    ('import time; time.sleep(10)', .05, 'timeout_or_deadline')])
def test_failure_and_timeout_never_mark_success(tmp_path, monkeypatch, code, timeout, status):
    import json
    monkeypatch.setattr(runner, 'commands', lambda *args: [
        ('july15', 'prepare', [sys.executable, '-c', code])])
    output = tmp_path / 'run'
    assert not runner.run(tmp_path, output,
        datetime.now(timezone.utc)+timedelta(hours=1), timeout_seconds=timeout)
    assert json.loads((output / 'stage_status.json').read_text())[0]['status'] == status
    assert (output / 'PLAN.json').exists()
    assert not (output / 'SUCCESS').exists()


def test_exit_zero_without_child_success_is_not_success(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, 'commands', lambda *args: [
        ('july15', 'prepare', [sys.executable, '-c', 'pass'])])
    output = tmp_path / 'run'
    with pytest.raises(AssertionError, match='Missing'):
        runner.run(tmp_path, output, datetime.now(timezone.utc)+timedelta(hours=1))
    assert not (output / 'SUCCESS').exists()
