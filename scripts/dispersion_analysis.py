#!/usr/bin/env python3
"""
Theoretical dispersion relation analysis for the 1D wave equation
discretized with finite elements and various time-stepping schemes.

This script computes and plots the ratio ω_h / ω_exact as a function of
the normalized wavenumber ξ = kh ∈ (0, π) for:

  - Central difference + lumped mass
  - Central difference + consistent mass
  - Newmark (β=1/4, γ=1/2) + lumped mass
  - Newmark (β=1/4, γ=1/2) + consistent mass

The 1D analysis gives insight into phase-speed errors (dispersion)
and is directly relevant to the 2D solver behaviour.

Mass and stiffness entries for P1 elements on a uniform 1D mesh of
spacing h:

  Consistent mass:  M_ij = h/6 * [1, 4, 1]   (tridiagonal Toeplitz)
  Lumped mass:       M_L  = h * I
  Stiffness:         K_ij = 1/h * [-1, 2, -1]

For a Fourier mode u_j = exp(i k j h), the eigenvalue of the generalized
problem K φ = λ M φ is:

  λ_K(ξ) = (2/h)(1 − cos ξ),    ξ = kh
  λ_M^L(ξ) = h
  λ_M^C(ξ) = (h/3)(2 + cos ξ)

The exact dispersion relation:  ω = k = ξ/h.
The spatial semi-discrete relation:  ω_h² = λ_K / λ_M.

For central difference:
  2(1 − cos(ω_h Δt)) / Δt² = λ_K / λ_M
  → cos(ω_h^d Δt) = 1 − (Δt²/2) * λ_K / λ_M
  → ω_h^d = (1/Δt) arccos(1 − (Δt²/2) ω_h²)

For Newmark (β=1/4, γ=1/2) — trapezoidal rule:
  The scheme is energy-conserving and the effective dispersion satisfies
  the same semi-discrete relation but solved implicitly.

Usage:
    python3 dispersion_analysis.py [result_dir] [--show]

If result_dir is given, saves to result_dir/dispersion_relation.png.
Otherwise saves to current directory.
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np


def main():
    args = [a for a in sys.argv[1:] if a != '--show']
    show_plot = '--show' in sys.argv

    out_dir = args[0] if args else '.'

    # Normalized wavenumber ξ = kh ∈ (0, π)
    xi = np.linspace(0.01, np.pi, 500)

    # ---- Spatial eigenvalues (P1 FEM, 1D uniform mesh) ----
    # Stiffness eigenvalue: λ_K = (2/h)(1 − cos ξ)   → normalized by h: 2(1−cos ξ)
    lam_K = 2.0 * (1.0 - np.cos(xi))

    # Lumped mass eigenvalue: λ_M = h  → normalized: 1
    lam_M_lumped = np.ones_like(xi)

    # Consistent mass eigenvalue: λ_M = (h/3)(2 + cos ξ) → normalized by h: (2+cos ξ)/3
    lam_M_consistent = (2.0 + np.cos(xi)) / 3.0

    # Semi-discrete ω²_h (spatial only)  for each mass type
    omega2_lumped = lam_K / lam_M_lumped       # = 2(1−cos ξ)
    omega2_consistent = lam_K / lam_M_consistent

    # ω_h * h  (semi-discrete, spatial only)
    omega_h_lumped = np.sqrt(omega2_lumped)
    omega_h_consistent = np.sqrt(omega2_consistent)

    # Exact: ω_exact * h = ξ
    omega_exact = xi

    # Ratio: ω_h / ω_exact = ω_h*h / ξ
    ratio_semi_lumped = omega_h_lumped / omega_exact
    ratio_semi_consistent = omega_h_consistent / omega_exact

    # ---- Fully-discrete dispersion (with time stepping) ----
    # CFL numbers to display
    cfl_values = [0.5, 0.8, 1.0]

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # --- Panel 1: Semi-discrete (spatial discretization only) ---
    ax = axes[0]
    ax.plot(xi / np.pi, ratio_semi_lumped, '-', color='#2196F3', linewidth=2,
            label='Lumped mass')
    ax.plot(xi / np.pi, ratio_semi_consistent, '-', color='#E91E63', linewidth=2,
            label='Consistent mass')
    ax.axhline(1.0, color='black', linewidth=0.8, linestyle='--', alpha=0.5, label='Exact')
    ax.set_xlabel('$\\xi / \\pi = kh/\\pi$', fontsize=13)
    ax.set_ylabel('$\\omega_h / \\omega_{\\rm exact}$', fontsize=13)
    ax.set_title('Semi-discrete dispersion\n(spatial only, P1 FEM)', fontsize=13)
    ax.set_xlim(0, 1)
    ax.set_ylim(0.5, 1.15)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11)

    # --- Panel 2: Central difference (fully discrete) ---
    ax = axes[1]
    for cfl_idx, cfl in enumerate(cfl_values):
        for mass_label, omega2_h in [('lumped', omega2_lumped),
                                      ('consistent', omega2_consistent)]:
            # σ = Δt * ω_h = cfl * h * ω_h / h  (but h cancels in normalized form)
            # Actually Δt = cfl * h, so Δt² ω²_h = cfl² * omega2_h (already normalized by h)
            arg = 1.0 - 0.5 * cfl**2 * omega2_h

            # Stability: |arg| ≤ 1
            valid = np.abs(arg) <= 1.0
            omega_d = np.full_like(xi, np.nan)
            omega_d[valid] = np.arccos(arg[valid]) / cfl  # normalized by h

            ratio = omega_d / omega_exact

            style = '-' if mass_label == 'lumped' else '--'
            color = ['#2196F3', '#4CAF50', '#FF9800'][cfl_idx]
            ax.plot(xi[valid] / np.pi, ratio[valid], linestyle=style, color=color,
                    linewidth=1.5, label=f'CFL={cfl}, {mass_label}')

    ax.axhline(1.0, color='black', linewidth=0.8, linestyle='--', alpha=0.5)
    ax.set_xlabel('$\\xi / \\pi = kh/\\pi$', fontsize=13)
    ax.set_ylabel('$\\omega_h^d / \\omega_{\\rm exact}$', fontsize=13)
    ax.set_title('Central difference\n(fully discrete, P1 FEM)', fontsize=13)
    ax.set_xlim(0, 1)
    ax.set_ylim(0.5, 1.15)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2)

    # --- Panel 3: Newmark β=1/4, γ=1/2 (fully discrete) ---
    ax = axes[2]
    beta = 0.25
    for cfl_idx, cfl in enumerate(cfl_values):
        for mass_label, omega2_h in [('lumped', omega2_lumped),
                                      ('consistent', omega2_consistent)]:
            # Newmark (trapezoidal):
            # (1 + β Δt² ω²_h) ω̃² = ω²_h
            # cos(ω̃ Δt) = (2 − (1−2β) Δt² ω²_h) / (2 + 2β Δt² ω²_h)
            # with β=1/4:
            # cos(ω̃ Δt) = (2 − 0.5 Δt² ω²_h) / (2 + 0.5 Δt² ω²_h)
            dt2_omega2 = cfl**2 * omega2_h
            numer = 2.0 - (1.0 - 2.0 * beta) * dt2_omega2
            denom = 2.0 + 2.0 * beta * dt2_omega2

            arg = numer / denom
            valid = np.abs(arg) <= 1.0
            omega_d = np.full_like(xi, np.nan)
            omega_d[valid] = np.arccos(arg[valid]) / cfl

            ratio = omega_d / omega_exact

            style = '-' if mass_label == 'lumped' else '--'
            color = ['#2196F3', '#4CAF50', '#FF9800'][cfl_idx]
            ax.plot(xi[valid] / np.pi, ratio[valid], linestyle=style, color=color,
                    linewidth=1.5, label=f'CFL={cfl}, {mass_label}')

    ax.axhline(1.0, color='black', linewidth=0.8, linestyle='--', alpha=0.5)
    ax.set_xlabel('$\\xi / \\pi = kh/\\pi$', fontsize=13)
    ax.set_ylabel('$\\omega_h^d / \\omega_{\\rm exact}$', fontsize=13)
    ax.set_title('Newmark ($\\beta=1/4, \\gamma=1/2$)\n(fully discrete, P1 FEM)', fontsize=13)
    ax.set_xlim(0, 1)
    ax.set_ylim(0.5, 1.15)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2)

    fig.suptitle('Dispersion relation: numerical vs exact phase velocity', fontsize=15, y=1.03)
    fig.tight_layout()

    out = os.path.join(out_dir, 'dispersion_relation.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Saved {out}")

    # ---- Additional plot: dissipation (amplitude ratio) ----
    fig2, axes2 = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes2[0]
    ax.axhline(1.0, color='black', linewidth=0.8, linestyle='--', alpha=0.5, label='Exact (no dissipation)')
    ax.set_xlabel('$\\xi / \\pi$', fontsize=13)
    ax.set_ylabel('Amplification factor $|G|$', fontsize=13)
    ax.set_title('Central difference\namplification per step', fontsize=13)

    for cfl_idx, cfl in enumerate(cfl_values):
        for mass_label, omega2_h in [('lumped', omega2_lumped),
                                      ('consistent', omega2_consistent)]:
            # For central difference, |G| = 1 exactly (energy conserving)
            # G = exp(±i ω_d Δt), so |G| = 1 when ω_d is real
            arg = 1.0 - 0.5 * cfl**2 * omega2_h
            valid = np.abs(arg) <= 1.0
            amp = np.ones_like(xi)
            amp[~valid] = np.nan

            style = '-' if mass_label == 'lumped' else '--'
            color = ['#2196F3', '#4CAF50', '#FF9800'][cfl_idx]
            ax.plot(xi[valid] / np.pi, amp[valid], linestyle=style, color=color,
                    linewidth=1.5, label=f'CFL={cfl}, {mass_label}')

    ax.set_xlim(0, 1)
    ax.set_ylim(0.95, 1.05)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    ax.text(0.5, 0.15, 'Central difference: $|G|=1$ exactly\n(no numerical dissipation)',
            transform=ax.transAxes, fontsize=11, ha='center',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    ax = axes2[1]
    ax.axhline(1.0, color='black', linewidth=0.8, linestyle='--', alpha=0.5, label='Exact')
    ax.set_xlabel('$\\xi / \\pi$', fontsize=13)
    ax.set_ylabel('Amplification factor $|G|$', fontsize=13)
    ax.set_title('Newmark ($\\beta=1/4, \\gamma=1/2$)\namplification per step', fontsize=13)

    for cfl_idx, cfl in enumerate(cfl_values):
        for mass_label, omega2_h in [('lumped', omega2_lumped),
                                      ('consistent', omega2_consistent)]:
            # Newmark β=1/4, γ=1/2 is also non-dissipative (|G|=1)
            amp = np.ones_like(xi)
            style = '-' if mass_label == 'lumped' else '--'
            color = ['#2196F3', '#4CAF50', '#FF9800'][cfl_idx]
            ax.plot(xi / np.pi, amp, linestyle=style, color=color,
                    linewidth=1.5, label=f'CFL={cfl}, {mass_label}')

    ax.set_xlim(0, 1)
    ax.set_ylim(0.95, 1.05)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    ax.text(0.5, 0.15, 'Newmark ($\\beta$=1/4, $\\gamma$=1/2): $|G|=1$\n(no numerical dissipation)',
            transform=ax.transAxes, fontsize=11, ha='center',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    fig2.suptitle('Numerical dissipation analysis: amplification factor per time step',
                  fontsize=14, y=1.02)
    fig2.tight_layout()

    out2 = os.path.join(out_dir, 'dissipation_amplification.png')
    fig2.savefig(out2, dpi=150, bbox_inches='tight')
    print(f"Saved {out2}")

    if show_plot:
        plt.show()
    plt.close('all')


if __name__ == '__main__':
    main()
