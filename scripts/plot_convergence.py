#!/usr/bin/env python3
"""
Convergence study plots (h-refinement and dt-refinement).

Scans all result directories to extract final L2/H1 errors, then produces
log-log convergence plots with fitted convergence rates.

Usage:
    python3 plot_convergence.py <result_dir> [--show]
"""

import csv
import glob
import math
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np

# ------------------------------------------------------------------ helpers

def parse_run_dir_name(name):
    """
    Extract configuration from directory name like:
      mesh-mesh_square_h0.05-mode-eigenmode-time-cd-mass-lumped-p1-dt-0.005-T-2-...
    Returns dict with keys: h, dt, scheme, mass, p  (strings / floats).
    """
    info = {}
    # h
    m = re.search(r'h([\d.]+)', name)
    if m:
        info['h'] = float(m.group(1))
    # dt
    m = re.search(r'-dt-([\d.eE+-]+)', name)
    if m:
        info['dt'] = float(m.group(1))
    # time scheme
    m = re.search(r'-time-(\w+)', name)
    if m:
        info['scheme'] = m.group(1)
    # mass
    m = re.search(r'-mass-(\w+)', name)
    if m:
        info['mass'] = m.group(1)
    # p (fe degree)
    m = re.search(r'-p(\d+)-', name)
    if m:
        info['p'] = int(m.group(1))
    return info


def read_final_errors(csv_path):
    """Read last row of error CSV and return (L2, H1)."""
    rows = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    if not rows:
        return None, None
    last = rows[-1]
    l2 = float(last['L2_error'])
    h1 = float(last.get('H1_error', 0))
    return l2, h1


def fit_slope(x, y):
    """Fit slope in log-log space."""
    lx = np.log(np.array(x))
    ly = np.log(np.array(y))
    A = np.vstack([lx, np.ones_like(lx)]).T
    slope, _ = np.linalg.lstsq(A, ly, rcond=None)[0]
    return slope

# ------------------------------------------------------------------ main

def main():
    args = [a for a in sys.argv[1:] if a != '--show']
    show_plot = '--show' in sys.argv

    if not args:
        print("Usage: python3 plot_convergence.py <result_dir> [--show]")
        sys.exit(1)

    result_dir = args[0]

    # Collect all result sub-directories
    subdirs = sorted(glob.glob(os.path.join(result_dir, '*')))
    subdirs = [d for d in subdirs if os.path.isdir(d)]

    # ===================== h-convergence =====================
    # Group by (scheme, mass, p, dt) → dict of h → (L2, H1)
    h_groups = {}
    for d in subdirs:
        info = parse_run_dir_name(os.path.basename(d))
        if not all(k in info for k in ('h', 'dt', 'scheme', 'mass', 'p')):
            continue
        # Find error CSV
        csvs = glob.glob(os.path.join(d, 'error-*.csv'))
        if not csvs:
            continue
        l2, h1 = read_final_errors(csvs[0])
        if l2 is None:
            continue
        key = (info['scheme'], info['mass'], info['p'], info['dt'])
        h_groups.setdefault(key, []).append((info['h'], l2, h1))

    # For h-convergence, pick groups where dt is the smallest (most refined in time)
    # so temporal error doesn't dominate
    h_conv_data = {}
    for key, entries in h_groups.items():
        scheme, mass, p, dt = key
        label_key = (scheme, mass, p)
        if label_key not in h_conv_data or dt < h_conv_data[label_key][0]:
            h_conv_data[label_key] = (dt, sorted(entries, key=lambda x: x[0]))

    if h_conv_data:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        for (scheme, mass, p), (dt, entries) in sorted(h_conv_data.items()):
            if len(entries) < 2:
                continue
            hs = [e[0] for e in entries]
            l2s = [e[1] for e in entries]
            h1s = [e[2] for e in entries]
            label = f"{scheme}-{mass}-P{p} (dt={dt})"
            slope_l2 = fit_slope(hs, l2s)
            slope_h1 = fit_slope(hs, h1s) if all(v > 0 for v in h1s) else 0

            axes[0].loglog(hs, l2s, 'o-', label=f"{label}  slope={slope_l2:.2f}")
            if all(v > 0 for v in h1s):
                axes[1].loglog(hs, h1s, 's-', label=f"{label}  slope={slope_h1:.2f}")

        # Reference slopes
        h_ref = np.array([0.025, 0.2])
        for ax, (expected_p1, expected_p2, title) in zip(
            axes,
            [(2, 3, '$L^2$ error vs $h$'), (1, 2, '$H^1$ error vs $h$')]
        ):
            y_ref = h_ref ** expected_p1 * 0.5
            ax.loglog(h_ref, y_ref, 'k--', alpha=0.4, label=f'$O(h^{expected_p1})$ ref')
            y_ref2 = h_ref ** expected_p2 * 0.3
            ax.loglog(h_ref, y_ref2, 'k:', alpha=0.4, label=f'$O(h^{expected_p2})$ ref')
            ax.set_xlabel('Mesh size $h$', fontsize=13)
            ax.set_ylabel('Error', fontsize=13)
            ax.set_title(title, fontsize=14)
            ax.grid(True, which='both', alpha=0.3)
            ax.legend(fontsize=9)
            ax.invert_xaxis()

        fig.suptitle('Spatial convergence (h-refinement)', fontsize=15, y=1.02)
        fig.tight_layout()
        out = os.path.join(result_dir, 'convergence_h.png')
        fig.savefig(out, dpi=150, bbox_inches='tight')
        print(f"Saved {out}")
        if show_plot:
            plt.show()
        plt.close(fig)

    # ===================== dt-convergence =====================
    dt_groups = {}
    for key, entries in h_groups.items():
        scheme, mass, p, dt = key
        # Group by (scheme, mass, p, h) → dict of dt → (L2, H1)
        for h_val, l2, h1 in entries:
            gk = (scheme, mass, p, h_val)
            dt_groups.setdefault(gk, []).append((dt, l2, h1))

    # Pick groups with smallest h → most refined in space
    dt_conv_data = {}
    for (scheme, mass, p, h_val), dts in dt_groups.items():
        label_key = (scheme, mass, p)
        if label_key not in dt_conv_data or h_val < dt_conv_data[label_key][0]:
            dt_conv_data[label_key] = (h_val, sorted(dts, key=lambda x: x[0]))

    if dt_conv_data:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        for (scheme, mass, p), (h_val, entries) in sorted(dt_conv_data.items()):
            if len(entries) < 2:
                continue
            dts_v = [e[0] for e in entries]
            l2s = [e[1] for e in entries]
            h1s = [e[2] for e in entries]
            label = f"{scheme}-{mass}-P{p} (h={h_val})"
            slope_l2 = fit_slope(dts_v, l2s) if len(dts_v) >= 2 else 0
            slope_h1 = fit_slope(dts_v, h1s) if len(dts_v) >= 2 and all(v > 0 for v in h1s) else 0

            axes[0].loglog(dts_v, l2s, 'o-', label=f"{label}  slope={slope_l2:.2f}")
            if all(v > 0 for v in h1s):
                axes[1].loglog(dts_v, h1s, 's-', label=f"{label}  slope={slope_h1:.2f}")

        dt_ref = np.array([0.002, 0.04])
        for ax, title in zip(axes, ['$L^2$ error vs $\\Delta t$', '$H^1$ error vs $\\Delta t$']):
            y2 = dt_ref ** 2 * 2.0
            ax.loglog(dt_ref, y2, 'k--', alpha=0.4, label='$O(\\Delta t^2)$ ref')
            ax.set_xlabel('$\\Delta t$', fontsize=13)
            ax.set_ylabel('Error', fontsize=13)
            ax.set_title(title, fontsize=14)
            ax.grid(True, which='both', alpha=0.3)
            ax.legend(fontsize=9)

        fig.suptitle('Temporal convergence ($\\Delta t$-refinement)', fontsize=15, y=1.02)
        fig.tight_layout()
        out = os.path.join(result_dir, 'convergence_dt.png')
        fig.savefig(out, dpi=150, bbox_inches='tight')
        print(f"Saved {out}")
        if show_plot:
            plt.show()
        plt.close(fig)

    if not h_conv_data and not dt_conv_data:
        print("[info] Not enough data for convergence plots.")


if __name__ == '__main__':
    main()
