import io
from pathlib import Path
import sys
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from frame_progress import FrameProgress


def test_resume_eta_uses_only_newly_processed_frames():
    clock = [100.]
    output = io.StringIO()
    progress = FrameProgress(1000, initial_done=600, stream=output, clock=lambda: clock[0])
    progress.update(600, phase='reading first resumed batch')
    assert '600/1000' in output.getvalue() and 'ETA --' in output.getvalue()
    clock[0] = 110.
    progress.update(620)
    # 20 newly processed frames / 10 s; 380 remaining -> 190 s.
    assert '62.0%' in output.getvalue() and 'ETA ~190s' in output.getvalue()
    progress.update(1000)
    assert output.getvalue().splitlines()[-1].count('ETA --') == 1


def test_new_run_eta_is_unchanged():
    clock = [0.]
    output = io.StringIO()
    progress = FrameProgress(100, stream=output, clock=lambda: clock[0])
    clock[0] = 5.
    progress.update(10)
    assert 'ETA ~45s' in output.getvalue()


@pytest.mark.parametrize('initial', [-1, 11])
def test_invalid_resume_offset_is_rejected(initial):
    with pytest.raises(ValueError, match='initial_done'):
        FrameProgress(10, initial_done=initial)
