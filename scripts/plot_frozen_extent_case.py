"""Inspect an extreme immutable seed/extent pair without local resegmentation."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import ndimage

from compare_frozen_extents import digest_frame
from compare_rain_thresholds import sha_file
from run_npy_validation import disk
from step.envelope import expand_seeded_envelope
from step.identification import identify_frame


def run(extents, output, metric='area_ratio'):
    if output.exists():
        raise FileExistsError(output)
    if not (extents / 'SUCCESS').is_file():
        raise ValueError('A successful fixed-seed extent run is required')
    config = json.loads((extents / 'configuration.json').read_text())
    source, metadata = Path(config['prepared_input']), Path(config['metadata'])
    if sha_file(source) != config['input_sha256'] or sha_file(metadata) != config['metadata_sha256']:
        raise ValueError('Prepared input identity changed')
    case = json.loads((extents / 'extreme_objects.json').read_text())[metric][0]
    frame, label = case['frame'], case['core_label']
    rain = np.array(np.load(source, mmap_mode='r')[frame], copy=True)
    settings = config['preset']['identification']
    seeds = identify_frame(rain, disk(settings['bridge_radius_cells']),
        min_size=settings['min_size_cells'], threshold=1)
    hashes = json.loads((extents / 'seed_frame_hashes.json').read_text())
    if digest_frame(seeds) != hashes[frame]:
        raise AssertionError('Selected seed raster changed')
    # Segment the FULL grid first; cropping before watershed would change ownership.
    high = expand_seeded_envelope(rain, seeds, .5)
    low = expand_seeded_envelope(rain, seeds, .1)
    selected = low if case['threshold_mm_h'] == .1 else high
    if case['threshold_mm_h'] not in (.1, .5) or int((selected == label).sum()) != case['envelope_cells']:
        raise AssertionError('Extreme case does not reproduce saved extent')
    yy, xx = np.nonzero(low == label)
    y0, y1 = max(0, int(yy.min())-12), min(rain.shape[0], int(yy.max())+13)
    x0, x1 = max(0, int(xx.min())-12), min(rain.shape[1], int(xx.max())+13)
    view = np.s_[y0:y1, x0:x1]
    cy, cx = np.nonzero(seeds == label)
    component, _ = ndimage.label(np.isfinite(rain) & (rain >= .1), np.ones((3, 3)))
    connected = np.unique(component[seeds == label])
    neighbors = np.unique(seeds[np.isin(component, connected)])
    details = dict(case, timestamp=json.loads(metadata.read_text())['timestamps'][frame],
        core_max_rain_mm_h=float(rain[seeds == label].max()),
        core_mean_rain_mm_h=float(rain[seeds == label].mean()),
        extent_mean_rain_mm_h=float(rain[low == label].mean()),
        extent_mean_rain_threshold_mm_h=.1,
        extent_0p5_cells=int((high == label).sum()),
        extent_0p1_cells=int((low == label).sum()),
        seed_labels_in_touched_weak_components=neighbors[neighbors > 0].tolist(),
        crop_yx=[y0, y1, x0, x1],
        note='Full-grid segmentation before plotting crop; diagnostic, not a false-assignment verdict')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    bounds = (x0-.5, x1-.5, y0-.5, y1-.5)
    image = axes[0].imshow(rain[view], origin='lower', extent=bounds,
                           cmap='Blues', vmin=0, vmax=3, interpolation='nearest')
    axes[0].set_title('Rain mm/h (display capped at 3)')
    fig.colorbar(image, ax=axes[0], shrink=.7, label='mm/h')
    palette = ListedColormap(['white', '#cce5ee', '#c0c0c0', '#f09a32'])
    for ax, partition, threshold in zip(axes[1:], (high, low), (.5, .1)):
        layer = np.zeros(rain.shape, dtype=np.uint8)
        layer[np.isfinite(rain) & (rain >= threshold)] = 1  # unassigned wet
        layer[partition > 0] = 2  # other seed owners
        layer[partition == label] = 3
        ax.imshow(layer[view], origin='lower', extent=bounds, cmap=palette,
                  vmin=0, vmax=3, interpolation='nearest')
        ax.set_title(f'{threshold} extent: orange = seed {label}')
    for ax in axes:
        ax.scatter([cx.mean()], [cy.mean()], c='red', marker='x', s=65,
                   label=f'seed {label} center', zorder=5)
        ax.set(xlabel='x (full-grid cells)', ylabel='y (full-grid cells)')
        ax.set_xlim(x0-.5, x1-.5)
        ax.set_ylim(y0-.5, y1-.5)
    axes[0].legend(loc='lower left', fontsize=8)
    fig.suptitle(f"{details['timestamp']} | node {case['core_node_id']}, branch {case['branch_id']} | "
        f"{case['core_cells']} seed cells -> {case['envelope_cells']} extent cells at {case['threshold_mm_h']} mm/h\n"
        'gray: other seed ownership; pale blue: unassigned wet; red cross: selected seed center')
    output.mkdir(parents=True)
    fig.savefig(output / 'case.png', dpi=150)
    plt.close(fig)
    (output / 'case.json').write_text(json.dumps(details, indent=2)+'\n')
    return details


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('extents', type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--metric', choices=('area_ratio', 'centroid_drift', 'any_seed_distance'),
                        default='area_ratio')
    args = parser.parse_args()
    print(json.dumps(run(args.extents, args.output_dir, args.metric), indent=2))
