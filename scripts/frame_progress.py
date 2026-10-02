"""Dependency-free frame progress readable in interactive and batch logs."""
import sys
import time


class FrameProgress:
    def __init__(self, total, enabled=True, stream=None, clock=time.monotonic,
                 initial_done=0):
        if total <= 0:
            raise ValueError('total must be positive')
        if not 0 <= initial_done <= total:
            raise ValueError('initial_done must be between zero and total')
        self.total, self.done = total, initial_done
        self.initial_done = initial_done
        self.enabled = enabled
        self.stream = stream if stream is not None else sys.stderr
        self.clock, self.started = clock, clock()
        self.unit = 'frames'
        self.verb = 'tracked'

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
        processed = self.done - self.initial_done
        eta = (f'~{elapsed/processed*(self.total-self.done):.0f}s'
               if processed > 0 and self.done < self.total else '--')
        print(f'[{bar}] {self.done}/{self.total} {self.unit} {self.verb} '
              f'({fraction:.1%}) | elapsed {elapsed:.0f}s | ETA {eta} | {phase}',
              file=self.stream, flush=True)

    def finish(self, success):
        self.update(phase='SUCCESS: outputs saved' if success else
                    'FAILED: see traceback; outputs may be partial')
