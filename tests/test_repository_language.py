"""Prevent Chinese prose from returning to the English-only beta snapshot."""

from pathlib import Path
import re


def test_repository_text_has_no_han_characters():
    root = Path(__file__).resolve().parents[1]
    text_suffixes = {'.md', '.txt', '.py', '.toml', '.json', '.slurm'}
    paths = [p for p in root.iterdir() if p.is_file() and p.suffix in text_suffixes]
    for directory in ('docs', 'configs', 'scripts', 'step', 'tests'):
        paths.extend(p for p in (root / directory).rglob('*')
                     if p.is_file() and p.suffix in text_suffixes)
    han = re.compile('[\u3400-\u4dbf\u4e00-\u9fff\U00020000-\U000323af]')
    offenders = []
    for path in sorted(paths):
        for number, line in enumerate(path.read_text(encoding='utf-8').splitlines(), 1):
            if han.search(line):
                offenders.append(f'{path.relative_to(root)}:{number}')
    assert not offenders, 'English-only repository policy violated: ' + ', '.join(offenders)
