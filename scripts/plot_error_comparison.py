#!/usr/bin/env python3
"""
Error comparison (dispersion analysis via error growth).

Overlays L2 / H1 error time-histories from multiple runs to show how
numerical dispersion causes phase error accumulation over time.

Usage:
    python3 plot_error_comparison.py <result_dir> [--show]
"""

import csv
import glob
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np

COLORS = ['#2196F3', '#E91E63', '#4CAF50', '#FF9800', '#9C27B0', '#00BCD4',
          '#795548', '#607D8B']
STYLES = ['-', '--', '-.', ':']


def parse_run_dir_name(name):
    info = {}
    m = re.search(r'h([\d.]+)', name)
    if m: info['h'] = float(m.group(1))
    m = re.search(r'-dt-([\d.eE+-]+)', name)
    if m: info['dt'] = float(m.group(1))
    m = re.search(r'-time-(\w+)', name)
    if m: info['scheme'] = m.group(1)
    m = re.search(r'-mass-(\w+)', name)
    if m: info['mass'] = m.group(1)
    m = re.search(r'-p(\d+)-', name)
    if m: info['p'] = int(m.group(1))
    return info


def read_error_csv(csv_path):
    times, l2s, h1s = [], [], []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            times.append(float(row['time']))
            l2s.append(float(row['L2_error']))
            h1v = row.get('H1_error', '')
            h1s.append(float(h1v) if h1v else 0.0)
    return np.array(times), np.array(l2s), np.array(h1s)


def make_label(info):
    s = info.get('scheme', '?')
    m = info.get('mass', '?')
    p = info.get('p', 1)
    return f"{s}-{m}-P{p}"


def main():
    args = [a for a in sys.argv[1:] if a != '--show']
    show_plot = '--show' in sys.argv

    if not args:
        print("Usage: python3 plot_error_comparison.py <result_dir> [--show]")
        sys.exit(1)

    result_dir = args[0]
    subdirs = sorted(glob.glob(os.path.join(result_dir, '*')))
    subdirs = [d for d in subdirs if os.path.isdir(d)]

    # Group by (h, dt) for scheme comparison
    groups = {}
    for d in subdirs:
        info = parse_run_dir_name(os.path.basename(d))
        if not all(k in info for k in ('h', 'dt', 'scheme', 'mass')):
            continue
        csvs = glob.glob(os.path.join(d, 'error-*.csv'))
        if not csvs:
            continue
        key = (info['h'], info['dt'])
        groups.setdefault(key, []).append((info, csvs[0]))

    comparison_groups = {k: v for k, v in groups.items() if len(v) >= 2}

    if not comparison_groups:
        print("[info] No scheme comparison data found.")
        return

    for (h, dt), entries in sorted(comparison_groups.items()):
        fig, axes = plt.subplots(2, 2, figsize=(15, 11))

        for idx, (info, csv_path) in enumerate(entries):
            times, l2s, h1s = read_error_csv(csv_path)
            label = make_label(info)
            color = COLORS[idx % len(COLORS)]
            style = STYLES[idx % len(STYLES)]

            # L2 error (linear)
            axes[0, 0].plot(times, l2s, color=color, linestyle=style,
                            linewidth=1.5, label=label)
            # L2 error (log)
            axes[0, 1].semilogy(times, l2s, color=color, linestyle=style,
                                linewidth=1.5, label=label)
            # H1 error (linear)
            if np.any(h1s > 0):
                axes[1, 0].plot(times, h1s, color=color, linestyle=style,
                                linewidth=1.5, label=label)
                # H1 error (log)
                axes[1, 1].semilogy(times, h1s, color=color, linestyle=style,
                                    linewidth=1.5, label=label)

        titles = [
            ('$L^2$ error vs time', 'Error'),
            ('$L^2$ error vs time (log)', 'Error (log)'),
            ('$H^1$ error vs time', 'Error'),
            ('$H^1$ error vs time (log)', 'Error (log)'),
        ]
        for ax, (title, ylabel) in zip(axes.flat, titles):
            ax.set_xlabel('Time $t$', fontsize=12)
            ax.set_ylabel(ylabel, fontsize=12)
            ax.set_title(title, fontsize=13)
            ax.grid(True, which='both', alpha=0.3)
            ax.legend(fontsize=9)

        fig.suptitle(
            f'Error comparison — Numerical dispersion  '
            f'($h={h}$, $\\Delta t={dt}$)',
            fontsize=14, y=1.02)
        fig.tight_layout()

        out = os.path.join(result_dir, f'error_comparison_h{h}_dt{dt}.png')
        fig.savefig(out, dpi=150, bbox_inches='tight')
        print(f"Saved {out}")
        if show_plot:
            plt.show()
        plt.close(fig)


if __name__ == '__main__':
    main()
