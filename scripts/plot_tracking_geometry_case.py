"""Review an endpoint change; masks are segmented on the full grid before cropping."""
import argparse
import json
from pathlib import Path

import numpy as np

from audit_reference_stream import inspect
from compare_frozen_extents import digest_frame, manifest
from compare_rain_thresholds import dump, sha_file
from run_npy_validation import disk
from step.envelope import expand_seeded_envelope
from step.identification import identify_frame
from step.tracking import _objects


def select_case(seed_objects, changes):
    seed = {(int(r['time']), int(r['local_label'])): r for r in seed_objects}
    if not changes:
        raise ValueError('No changes in requested class')
    def rank(r):
        endpoints = [seed[tuple(r['pair'][i:i+2])] for i in (0, 2)]
        return (min(float(v['precipitation_sum']) for v in endpoints), tuple(r['pair']))
    return max(changes, key=rank)


def run(experiment, output, geometry='extent_0p1', kind='added'):
    if output.exists():
        raise FileExistsError(output)
    if not (experiment / 'SUCCESS').exists():
        raise ValueError('Successful experiment required')
    plan = json.loads((experiment / 'PLAN.json').read_text())
    extents = Path(plan['source'])
    config = json.loads((extents / 'configuration.json').read_text())
    if manifest(Path(config['reference'])) != plan['source_graph_files_sha256']:
        raise ValueError('Original graph identity changed')
    source, metadata = Path(config['prepared_input']), Path(config['metadata'])
    if sha_file(source) != config['input_sha256'] or sha_file(metadata) != config['metadata_sha256']:
        raise ValueError('Prepared input identity changed')
    _, objects, _, _, _ = inspect(Path(config['reference']))
    changes = json.loads((experiment / (geometry+'_edge_changes.json')).read_text())[kind]
    case = select_case(objects, changes)
    by_key = {(int(r['time']), int(r['local_label'])): int(r['node_id']) for r in objects}
    parent_id, child_id = (by_key[tuple(case['pair'][i:i+2])] for i in (0, 2))
    context = {}
    for name in ('seed_control', 'extent_0p5', 'extent_0p1'):
        _, _, edges, _, _ = inspect(experiment / name)
        context[name] = [r for r in edges if int(r['parent_node_id']) == parent_id or
                        int(r['child_node_id']) == child_id]
    hashes = json.loads((extents / 'seed_frame_hashes.json').read_text())
    rain = np.load(source, mmap_mode='r')
    settings = config['preset']['identification']
    frames, layers, stats = [], [], []
    for frame, label in (case['pair'][:2], case['pair'][2:]):
        values = np.array(rain[frame], copy=True)
        seeds = identify_frame(values, disk(settings['bridge_radius_cells']),
            min_size=settings['min_size_cells'], threshold=1)
        if digest_frame(seeds) != hashes[frame]:
            raise AssertionError('Saved seed identity changed')
        partitions = [seeds, expand_seeded_envelope(values, seeds, .5),
                      expand_seeded_envelope(values, seeds, .1)]
        masks = [part == label for part in partitions]
        measured = [{obj.label: obj for obj in _objects(part, values)}[label] for part in partitions]
        frames.append(values)
        layers.append(masks)
        stats.append([dict(area=v.area, centroid=list(v.centroid), mean=v.mean_intensity,
                           max=v.max_intensity) for v in measured])
    yy, xx = np.nonzero(np.logical_or.reduce([m for row in layers for m in row]))
    ny, nx = rain.shape[1:]
    y0, y1 = max(0, int(yy.min())-12), min(ny, int(yy.max())+13)
    x0, x1 = max(0, int(xx.min())-12), min(nx, int(xx.max())+13)
    view, bounds = np.s_[y0:y1, x0:x1], (x0-.5, x1-.5, y0-.5, y1-.5)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    fig, axes = plt.subplots(2, 4, figsize=(15, 8), layout='constrained')
    stamps = json.loads(metadata.read_text())['timestamps']
    for i, (values, masks, measures) in enumerate(zip(frames, layers, stats)):
        frame, label = case['pair'][2*i:2*i+2]
        im = axes[i, 0].imshow(values[view], origin='lower', extent=bounds,
                             cmap='Blues', vmin=0, vmax=10)
        axes[i, 0].set_title(f'{stamps[frame]}\nRain mm/h; display cap 10')
        fig.colorbar(im, ax=axes[i, 0], shrink=.7)
        for ax, mask, measure, threshold in zip(axes[i, 1:], masks, measures, (1., .5, .1)):
            ax.imshow((values[view] >= threshold).astype(int) + 2*mask[view], origin='lower',
                extent=bounds, cmap=ListedColormap(['white', '#e2e2e2', '#f09a32', '#f09a32']),
                vmin=0, vmax=3, interpolation='nearest')
            cy, cx = measure['centroid']
            ax.scatter([cx], [cy], marker='x', color='red', s=40)
            ax.set_title(f"{threshold:g} geometry, seed {label}\n"
                         f"{measure['area']} cells, mean {measure['mean']:.2f}")
        for ax in axes[i]:
            ax.set(xlim=(x0-.5, x1-.5), ylim=(y0-.5, y1-.5),
                   xlabel='x full-grid cells', ylabel='y full-grid cells')
    fig.suptitle(f"{geometry}: {kind} endpoint pair {case['pair']}\n"
        'Selected by largest minimum original seed rain; not a representative or accuracy sample.\n'
        'Orange: selected owner; gray: other wet pixels; red: rain-weighted feature centroid.')
    output.mkdir(parents=True)
    fig.savefig(output / 'case.png', dpi=120)
    plt.close(fig)
    dump(output / 'case.json', dict(case=case, geometry=geometry, kind=kind,
        parent_node_id=parent_id, child_node_id=child_id, incident_edge_context=context,
        feature_statistics=stats, crop_yx=[y0, y1, x0, x1],
        selection='Largest minimum original seed precipitation_sum, tie by immutable endpoint pair',
        note='Full-grid replay before plotting crop. Added/removed edge is not a truth label. '
             'Incident endpoint edges are included to expose alternative links and events; '
             'a removed pair need not mean an entirely broken storm trajectory.'))
    return case


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('experiment', type=Path)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--geometry', choices=('extent_0p5', 'extent_0p1'), default='extent_0p1')
    p.add_argument('--kind', choices=('added', 'removed'), default='added')
    a = p.parse_args()
    print(run(a.experiment, a.output_dir, a.geometry, a.kind))
