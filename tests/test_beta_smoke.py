import json
from pathlib import Path
import sys

import pytest
pytest.importorskip('xarray')
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import beta_smoke


def test_beta_check_crosses_month_and_replays_gap_identity(tmp_path, monkeypatch):
    monkeypatch.setattr(beta_smoke, 'version', lambda name: beta_smoke.BETA_VERSION)
    out = tmp_path / 'smoke'
    result = beta_smoke.run(out, workers=2)
    assert all(result['checks'].values())
    assert result['frames'] == 12 and result['nodes'] == 11 and result['gap_edges'] == 1
    contract = json.loads((out / 'control/execution_contract.json').read_text())
    assert contract['snapshots'][0]['timestamp'].startswith('1996-01-31')
    assert contract['snapshots'][-1]['timestamp'].startswith('1996-02-01')
    assert (out / 'SUCCESS').exists()
    with pytest.raises(FileExistsError):
        beta_smoke.run(out)


def test_beta_check_rejects_stale_install_and_invalid_workers(tmp_path, monkeypatch):
    monkeypatch.setattr(beta_smoke, 'version', lambda name: '0.3.0')
    with pytest.raises(ValueError, match='reinstall'):
        beta_smoke.run(tmp_path / 'stale')
    assert not (tmp_path / 'stale').exists()
    with pytest.raises(ValueError, match='workers'):
        beta_smoke.run(tmp_path / 'workers', workers=0)
    assert not (tmp_path / 'workers').exists()
