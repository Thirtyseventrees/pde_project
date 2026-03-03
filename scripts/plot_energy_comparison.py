#!/usr/bin/env python3
"""
Energy comparison (dissipation analysis).

Overlays energy time-histories and relative energy drift from multiple runs
on a single figure to compare numerical dissipation across schemes.

Usage:
    python3 plot_energy_comparison.py <result_dir> [--show]

It looks for the "Experiment 3" runs (same mesh/dt, different scheme+mass)
by matching directory names.
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


def read_energy_csv(csv_path):
    times, energies = [], []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            times.append(float(row['time']))
            energies.append(float(row['energy']))
    return np.array(times), np.array(energies)


def make_label(info):
    s = info.get('scheme', '?')
    m = info.get('mass', '?')
    p = info.get('p', 1)
    return f"{s}-{m}-P{p}"


def main():
    args = [a for a in sys.argv[1:] if a != '--show']
    show_plot = '--show' in sys.argv

    if not args:
        print("Usage: python3 plot_energy_comparison.py <result_dir> [--show]")
        sys.exit(1)

    result_dir = args[0]
    subdirs = sorted(glob.glob(os.path.join(result_dir, '*')))
    subdirs = [d for d in subdirs if os.path.isdir(d)]

    # Group by (h, dt) → list of (info, energy_csv)
    groups = {}
    for d in subdirs:
        info = parse_run_dir_name(os.path.basename(d))
        if not all(k in info for k in ('h', 'dt', 'scheme', 'mass')):
            continue
        csvs = glob.glob(os.path.join(d, 'energy-*.csv'))
        if not csvs:
            continue
        key = (info['h'], info['dt'])
        groups.setdefault(key, []).append((info, csvs[0]))

    # Find groups with multiple entries (scheme comparison)
    comparison_groups = {k: v for k, v in groups.items() if len(v) >= 2}

    if not comparison_groups:
        print("[info] No scheme comparison data found (need multiple runs with same h,dt).")
        return

    for (h, dt), entries in sorted(comparison_groups.items()):
        fig, axes = plt.subplots(1, 2, figsize=(15, 6))

        for idx, (info, csv_path) in enumerate(entries):
            times, energies = read_energy_csv(csv_path)
            label = make_label(info)
            color = COLORS[idx % len(COLORS)]
            style = STYLES[idx % len(STYLES)]

            # Absolute energy
            axes[0].plot(times, energies, color=color, linestyle=style,
                         linewidth=1.5, label=label)

            # Relative drift
            E0 = energies[0] if energies[0] != 0 else 1e-16
            drift = (energies - E0) / E0 * 100.0
            axes[1].plot(times, drift, color=color, linestyle=style,
                         linewidth=1.5, label=label)

        axes[0].set_xlabel('Time $t$', fontsize=13)
        axes[0].set_ylabel('Discrete energy $E(t)$', fontsize=13)
        axes[0].set_title('Energy vs time', fontsize=14)
        axes[0].grid(True, alpha=0.3)
        axes[0].legend(fontsize=10)
        axes[0].ticklabel_format(useOffset=False)

        axes[1].set_xlabel('Time $t$', fontsize=13)
        axes[1].set_ylabel('$(E(t)-E_0)/E_0$ [%]', fontsize=13)
        axes[1].set_title('Relative energy drift (dissipation)', fontsize=14)
        axes[1].grid(True, alpha=0.3)
        axes[1].axhline(0, color='black', linewidth=0.5)
        axes[1].legend(fontsize=10)

        fig.suptitle(
            f'Energy comparison — Numerical dissipation  '
            f'($h={h}$, $\\Delta t={dt}$)',
            fontsize=14, y=1.02)
        fig.tight_layout()

        out = os.path.join(result_dir, f'energy_comparison_h{h}_dt{dt}.png')
        fig.savefig(out, dpi=150, bbox_inches='tight')
        print(f"Saved {out}")
        if show_plot:
            plt.show()
        plt.close(fig)


if __name__ == '__main__':
    main()
