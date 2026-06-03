#!/usr/bin/env python3
"""Validation: paper t + scaled G → can we reproduce paper mobility?"""
import sys; sys.path.insert(0,'.')
import numpy as np
from naph_transport.stage5_mobility_assemble import compute_direction_projections, compute_A0_alpha
from naph_transport.stage6_mobility_integrate import compute_mobility_vs_T
from naph_transport.io_utils import load_checkpoint
from naph_transport.constants import CM1_TO_RADS, REF_T_VALUES_EV

t_paper = np.array([REF_T_VALUES_EV[k] for k in ['a','b','c','ac','ab','abc']])
s5 = load_checkpoint('data/results/stage5_mobility_params.npz')
omega = s5['omega_cm1'] * CM1_TO_RADS
G_our = s5['G_total']
R_proj = compute_direction_projections()
A0_p = compute_A0_alpha(t_paper, R_proj)
T_range = np.arange(10, 310, 10)

# Paper G_eff for 3 modes: 0.40, 0.17, 0.10 (freqs 58.8, 82.0, 108.5 cm-1)
# Our G_total for closest modes:
print("Paper vs our G at similar frequencies:")
for pf, pG in [(58.8, 0.40), (82.0, 0.17), (108.5, 0.10)]:
    idx = np.argmin(np.abs(s5['omega_cm1'] - pf))
    print(f"  Paper: {pf:.0f} cm-1, G={pG:.2f}")
    print(f"  Our:   {s5['omega_cm1'][idx]:.0f} cm-1, G={G_our[idx]:.4f}")
    print(f"  Ratio G_paper/G_our: {pG/G_our[idx]:.1f}x")
    print()

# Test with scaled G
print(f'{"Scale":>6s} {"μ_c(300K)":>10s} {"Band":>8s} {"Hop":>8s} {"Mechanism":>18s} {"Paper match?":>15s}')
print('-'*75)

for scale in [1.0, 2.0, 3.0, 5.0, 10.0]:
    G_scaled = G_our * scale
    r = compute_mobility_vs_T(T_range, A0_p['c_prime'], s5['B_c_prime'], G_scaled, omega, verbose=False)
    mu = r['mu_total'][-1]
    band = r['mu_band'][-1]
    hop = r['mu_hopping'][-1]
    mech = 'Hopping' if hop > band else 'Band'
    match = '✓ YES!' if 0.05 < mu < 0.20 else ''
    print(f'{scale:6.1f} {mu:10.4f} {band:8.4f} {hop:8.4f} {mech:>18s} {match:>15s}')

# Best scale: ~3-5x gives μ ≈ 0.1-0.2
print(f'\nPaper G_eff/our G ratio: ~{(0.40+0.17+0.10)/(G_our[0]+G_our[1]+G_our[2]):.1f}x for 3 lowest modes')
print(f'Total ΣG_our: {np.sum(G_our):.2f}, ΣG_our×5: {np.sum(G_our)*5:.2f}')
