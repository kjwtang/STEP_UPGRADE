"""Exercise the documented prepared-rain API examples and their CSV schemas."""

import csv
from pathlib import Path
import re

import numpy as np


def test_readme_python_examples_and_outputs(tmp_path, monkeypatch):
    root = Path(__file__).resolve().parents[1]
    readme = (root / 'README.md').read_text()
    examples = re.findall(r'```python\n(.*?)\n```', readme, re.S)
    assert len(examples) == 6
    config_name = 'rainfall_lineage_4km_1h_reference_v1.json'
    (tmp_path / 'configs').mkdir()
    (tmp_path / 'configs' / config_name).write_text(
        (root / 'configs' / config_name).read_text())
    rain = np.zeros((4, 80, 100), dtype=np.float32)
    for frame in (0, 2, 3):
        rain[frame, 30:35, 20+frame:25+frame] = 2
    np.save(tmp_path / 'rain_mm_h.npy', rain)
    monkeypatch.chdir(tmp_path)
    monkeypatch.syspath_prepend(str(root / 'scripts'))
    namespace = {}
    for example in examples:
        if 'resume_state = load_tracking_state' in example:
            namespace['next_labels'] = namespace['labels'][-1:]
            namespace['next_rain'] = namespace['rain'][-1:]
        exec(compile(example, 'README.md', 'exec'), namespace)
    assert namespace['labels'].dtype == np.int32
    assert namespace['branches'].dtype == np.int64
    assert namespace['labels'].shape == rain.shape
    assert np.array_equal(namespace['frame_labels'], namespace['labels'][0])
    assert namespace['state'].last_time == 3
    assert namespace['next_state'].last_time == 4
    assert len(namespace['graph'].objects) == 3
    assert sum(edge.gap > 0 for edge in namespace['graph'].edges) == 1
    output = tmp_path / 'results/api_chunk_01'
    assert np.array_equal(np.load(output / 'tracked_labels.npy'), namespace['branches'])
    schemas = re.findall(r'^(?:time|event_id|branch_id),[^\n]+', readme, re.M)
    expected = {line.split(',')[1]: line.split(',') for line in schemas}
    for name, key in (('objects.csv', 'node_id'), ('edges.csv', 'parent_node_id'),
                      ('events.csv', 'time')):
        with (output / name).open(newline='') as handle:
            assert next(csv.reader(handle)) == expected[key]
    assert (output / 'state.json').is_file()
