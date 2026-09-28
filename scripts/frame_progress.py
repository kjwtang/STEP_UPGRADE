"""Dependency-free frame progress readable in interactive and batch logs."""
import sys
import time


class FrameProgress:
    def __init__(self, total, enabled=True, stream=None, clock=time.monotonic):
        self.total, self.done = total, 0
        self.enabled = enabled
        self.stream = stream if stream is not None else sys.stderr
        self.clock, self.started = clock, clock()
        self.unit = 'frames'

    def update(self, done=None, phase='tracking'):
        if done is not None:
            if not self.done <= done <= self.total:
                raise ValueError('Progress must be monotonic and bounded')
            self.done = done
        if not self.enabled:
            return
        elapsed = max(0., self.clock()-self.started)
        fraction = self.done/self.total
        filled = int(24*fraction)
        bar = '#' * filled + '-' * (24-filled)
        eta = (f'~{elapsed/self.done*(self.total-self.done):.0f}s'
               if 0 < self.done < self.total else '--')
        print(f'[{bar}] {self.done}/{self.total} {self.unit} tracked '
              f'({fraction:.1%}) | elapsed {elapsed:.0f}s | ETA {eta} | {phase}',
              file=self.stream, flush=True)

    def finish(self, success):
        self.update(phase='SUCCESS: outputs saved' if success else
                    'FAILED: see traceback; outputs may be partial')
