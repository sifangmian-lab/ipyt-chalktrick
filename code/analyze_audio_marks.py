"""Exploratory acoustic pulse timing and uncalibrated chalk-mark geometry.

Requires ffmpeg on PATH, numpy, scipy, matplotlib and Pillow. Original recordings
and photographs remain outside the repository. No sound-to-mark alignment is
assumed or fitted. Run with --help for input paths.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy import ndimage, signal


def write_csv(path, fields, rows):
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.writer(f)
        w.writerow(fields)
        w.writerows(rows)


def audio_envelope(a, band):
    fs = 48000
    b = signal.sosfiltfilt(signal.butter(4, band, fs=fs,
                                       btype='bandpass', output='sos'), a)
    env = np.sqrt(np.maximum(ndimage.gaussian_filter1d(b*b, fs*.00025), 0))[::12]
    # Local scaling preserves weak pulses without treating sound amplitude as force.
    local = ndimage.maximum_filter1d(env, 240)
    norm = env / (ndimage.gaussian_filter1d(local, 60) + 1e-4)
    return np.arange(len(env))/4000, env, norm


def detect(t, norm, window, min_distance, prominence):
    keep = (t >= window[0]) & (t <= window[1])
    ids = np.flatnonzero(keep)
    p, _ = signal.find_peaks(norm[keep], distance=round(min_distance*4000),
                             prominence=prominence, height=.20)
    return ids[p]


def summarize_intervals(times):
    dt = np.diff(times)
    return dict(n_pulses=len(times), median_interval_ms=float(np.median(dt)*1000),
                reciprocal_median_hz=float(1/np.median(dt)),
                interval_q10_q90_ms=(np.percentile(dt, [10, 90])*1000).tolist())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ['audio29', 'audio30', 'photo1', 'photo2', 'out']:
        parser.add_argument('--'+flag, type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    summary = {'date': '2026-10-29', 'pairing': {'29': 'photo1', '30': 'photo2'},
               'direction': 'left_to_right', 'pixel_coordinates': 'resized width=1824',
               'sound_to_mark_alignment': 'not established',
               'measurement_class': 'exploratory; freehand speed and load confounded',
               'sources': {}, 'audio': {}, 'photo': {}}
    for flag in ['audio29', 'audio30', 'photo1', 'photo2']:
        p = getattr(args, flag)
        summary['sources'][flag] = {'filename': p.name,
                                    'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}

    fig, axs = plt.subplots(2, 3, figsize=(16, 8), layout='constrained')
    colors = {'29': '#146e9e', '30': '#c05227'}
    settings = {'29': {'activity': (.22, 1.28), 'core': (.60, 1.08), 'distance': .012},
                '30': {'activity': (.63, .96), 'core': (.73, .92), 'distance': .006}}
    for row, label in enumerate(['29', '30']):
        cfg = settings[label]
        p = getattr(args, 'audio'+label)
        raw = subprocess.run(['ffmpeg', '-v', 'error', '-i', str(p), '-ac', '1',
                              '-ar', '48000', '-f', 'f32le', 'pipe:1'],
                             check=True, capture_output=True).stdout
        a = np.frombuffer(raw, dtype='<f4')
        t, env, norm = audio_envelope(a, (1200, 12000))
        events = detect(t, norm, cfg['activity'], cfg['distance'], .22)
        core = detect(t, norm, cfg['core'], cfg['distance'], .22)
        stats = summarize_intervals(t[core])
        stats.update(decoded_duration_s=len(a)/48000,
                     activity_window_s=cfg['activity'], core_window_s=cfg['core'],
                     activity_candidate_pulses=len(events),
                     effective_timing_grid_ms=.25,
                     confidence_interval='not estimated; quantiles describe within-run variation')
        sensitivity = []
        for band in [(600, 14000), (1200, 12000), (2000, 10000)]:
            tt, _, nn = audio_envelope(a, band)
            for prominence in [.18, .22, .28]:
                pp = detect(tt, nn, cfg['core'], cfg['distance'], prominence)
                sensitivity.append(dict(band_hz=band, prominence=prominence,
                                        **summarize_intervals(tt[pp])))
        stats['sensitivity'] = sensitivity
        summary['audio'][label] = stats
        times = t[events]
        write_csv(args.out/f'audio_{label}_candidate_events.csv',
                  ['event_index', 'time_s', 'interval_from_previous_ms', 'normalized_envelope',
                   'inside_core_window'],
                  [(i+1, ts, '' if i == 0 else (ts-times[i-1])*1000,
                    norm[ev], cfg['core'][0] <= ts <= cfg['core'][1])
                   for i, (ev, ts) in enumerate(zip(events, times))])
        write_csv(args.out/f'audio_{label}_core_events.csv', ['event_index', 'time_s'],
                  [(i+1, ts) for i, ts in enumerate(t[core])])
        # Keep a short common-duration excerpt to make the different rhythms visible.
        start = .80 if label == '29' else .75
        show = (t >= start) & (t <= start+.15)
        axs[row, 0].plot((t[show]-start)*1000, norm[show], color=colors[label], lw=1)
        axs[row, 0].set(title=f'{label} / photo {row+1}: acoustic pulse train',
                       xlabel='Elapsed time in excerpt (ms)', ylabel='Normalized envelope')
        mid = (times[1:]+times[:-1])/2
        axs[row, 1].plot(mid, np.diff(times)*1000, 'o-', ms=3,
                         color=colors[label], lw=1)
        axs[row, 1].axvspan(*cfg['core'], color=colors[label], alpha=.09)
        axs[row, 1].set(title=f'Core: {stats["median_interval_ms"]:.2f} ms '
                              f'(~{stats["reciprocal_median_hz"]:.0f} Hz)',
                       xlabel='Time in recording (s)', ylabel='Candidate pulse interval (ms)',
                       ylim=(0, 34))

    # Manual centerline guides, with background subtraction for geometry only.
    photo_cfg = {
        1: ([150, 300, 600, 900, 1200, 1500, 1800],
            [693, 678, 646, 619, 614, 605, 614], 4, (145, 1824), 24, 3, 2, 1805),
        2: ([0, 300, 600, 900, 1200, 1500, 1700],
            [606, 595, 604, 611, 609, 604, 607], 3, (0, 1720), 25, 6, 5, 1590)}
    for k in [1, 2]:
        p = getattr(args, 'photo'+str(k))
        with Image.open(p) as image:
            original_size = list(image.size)
            im = np.array(image.convert('RGB').resize((1824, 1368)))
        xg, yg, deg, xr, threshold, close, min_size, accepted_max = photo_cfg[k]
        x = np.arange(*xr)
        y = np.polyval(np.polyfit(xg, yg, deg), x)
        offsets = np.arange(-38, 39)
        coordinates = np.array([y[None, :]+offsets[:, None],
                                np.broadcast_to(x, (len(offsets), len(x)))])
        band = ndimage.map_coordinates(im.mean(axis=2), coordinates, order=1)
        residual = np.maximum(band-ndimage.gaussian_filter(band, (8, 16)), 0)
        profile = ndimage.gaussian_filter1d(np.mean(np.sort(residual, axis=0)[-15:], axis=0), .65)
        mask = ndimage.binary_closing(profile > threshold, structure=np.ones(close))
        labs, n = ndimage.label(mask)
        marks = []
        for lab in range(1, n+1):
            ids = np.flatnonzero(labs == lab)
            if len(ids) < min_size:
                continue
            weight = np.maximum(profile[ids]-threshold, 0)
            center = float(np.sum(x[ids]*weight)/weight.sum())
            status = 'interior_geometry' if center <= accepted_max else 'excluded_edge_or_overlap'
            marks.append((center, int(x[ids[0]]), int(x[ids[-1]]), status))
        write_csv(args.out/f'photo_{k}_candidate_geometry.csv',
                  ['candidate_index', 'center_x_px', 'left_px', 'right_px', 'status'],
                  [(i+1, *m) for i, m in enumerate(marks)])
        write_csv(args.out/f'photo_{k}_profile.csv',
                  ['x_px', 'local_brightness_excess_uncalibrated'], zip(x, profile))
        centers = np.array([m[0] for m in marks if m[3] == 'interior_geometry'])
        gap = np.diff(centers)
        axs[k-1, 2].plot((centers[1:]+centers[:-1])/2, gap, 'o-',
                         ms=3, lw=1, color=colors[str(28+k)])
        axs[k-1, 2].set(title=f'Photo {k}: spatial gaps (uncalibrated)',
                        xlabel='Image x position, left to right (px)', ylabel='Mark spacing (px)')
        summary['photo'][str(k)] = dict(original_size=original_size,
                                        candidate_segments=len(marks),
                                        interior_geometry_segments=len(centers),
                                        threshold=threshold,
                                        centerline_x=xg, centerline_y=yg,
                                        accepted_center_max_px=accepted_max,
                                        physical_scale='unknown; do not compare pixels across photos',
                                        total_mark_count='not certified; clipped/overlapping end marks')
        # A separate audit figure shows the segmentation against the actual band.
        fig2, ax2 = plt.subplots(2, 1, figsize=(15, 4), layout='constrained')
        ax2[0].imshow(band, cmap='gray', extent=(x[0], x[-1], 38, -38), aspect='auto')
        for i, m in enumerate(marks):
            ax2[0].text(m[0], -25, str(i+1), fontsize=7, color='red', ha='center')
        ax2[0].set(title=f'Photo {k}: candidate segments; end regions require manual review')
        ax2[1].plot(x, profile, lw=1)
        ax2[1].axhline(threshold, color='gray')
        ax2[1].set(xlabel='Resized image x coordinate (px)', ylabel='Local brightness excess')
        fig2.savefig(args.out/f'photo_{k}_segmentation_audit.png', dpi=140)
        plt.close(fig2)

    fig.suptitle('Exploratory paired recordings and chalk marks | 29 -> photo 1; 30 -> photo 2', fontsize=14)
    fig.savefig(args.out/'comparison.png', dpi=160)
    plt.close(fig)
    (args.out/'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    print(json.dumps({k: {field: summary['audio'][k][field]
                         for field in ['median_interval_ms', 'reciprocal_median_hz', 'interval_q10_q90_ms']}
                      for k in ['29', '30']}, ensure_ascii=False))


if __name__ == '__main__':
    main()
