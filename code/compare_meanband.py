#!/usr/bin/env python3
"""Compare mean-band vs 2x2 mobility with fixed units."""
import sys; sys.path.insert(0,'.')
import numpy as np
from naph_transport.stage5_mobility_assemble import (compute_direction_projections,
    compute_A0_alpha, compute_G_total, compute_B_alpha_lambda,
    dE_dQ_to_dimensionless_g)
from naph_transport.stage6_mobility_integrate import compute_mobility_vs_T
from naph_transport.io_utils import load_checkpoint, save_checkpoint
from naph_transport.constants import CM1_TO_RADS

t_mb = np.array([0.03427, -0.03133, -0.00537, -0.00782, 0.00473, 0.00069])
t_2x2 = np.array([0.03703, -0.02886, -0.00530, -0.00052, 0.06024, 0.01481])

s6_2x2 = load_checkpoint('data/results/stage6_mobility_results.npz')
s5 = load_checkpoint('data/results/stage5_mobility_params.npz')
s3 = load_checkpoint('data/results/stage3_eph_coupling.npz')
s2 = load_checkpoint('data/results/stage2_phonon.npz')

omega_cm1 = s2['selected_frequencies']
omega_rads = omega_cm1 * CM1_TO_RADS
d_eps = s3['d_epsilon0_dQ']; dt = s3['dt_dQ']
n_modes = len(omega_cm1)

g_loc = np.array([dE_dQ_to_dimensionless_g(d_eps[i],abs(omega_cm1[i])) for i in range(n_modes)])
g_nl = np.array([[dE_dQ_to_dimensionless_g(dt[i,j],abs(omega_cm1[i])) for j in range(6)] for i in range(n_modes)])

R = compute_direction_projections()
A_mb = compute_A0_alpha(t_mb, R)
A_2x2 = compute_A0_alpha(t_2x2, R)
G_tot = compute_G_total(g_loc, g_nl)
B = compute_B_alpha_lambda(g_nl, omega_cm1, R)
Tr = s6_2x2['T_range']

print("=" * 65)
print("  MEAN-BAND vs 2x2 MOBILITY (fixed units)")
print("=" * 65)
print("A0 ratios: ",end="")
for d in ['a','b','c_prime']: print(f"{d}={A_mb[d]/A_2x2[d]:.2f} ",end="")
print()

r_mb = {}
for d in ['a','b','c_prime']:
    r = compute_mobility_vs_T(Tr, A_mb[d], B[d], G_tot, omega_rads, verbose=False)
    for k,v in r.items(): r_mb[f'{k}_{d}'] = v

print(f"\n{'T':>5s} {'MB-a':>8s} {'2x2-a':>8s}  {'MB-b':>8s} {'2x2-b':>8s}  {'MB-c':>8s} {'2x2-c':>8s}")
print('-'*65)
for i,Tv in enumerate(Tr):
    if Tv in [10,50,100,200,300]:
        print(f'{Tv:5.0f} {r_mb["mu_total_a"][i]:8.2f} {s6_2x2["mu_total_a"][i]:8.2f}  '
              f'{r_mb["mu_total_b"][i]:8.2f} {s6_2x2["mu_total_b"][i]:8.2f}  '
              f'{r_mb["mu_total_c_prime"][i]:8.2f} {s6_2x2["mu_total_c_prime"][i]:8.2f}')

# Crossover
mb_b = r_mb['mu_band_c_prime']; mb_h = r_mb['mu_hopping_c_prime']
print(f"\nMean-band c' crossover:")
for Tv in [10,50,100,150,200,300]:
    i=np.argmin(np.abs(Tr-Tv))
    d='BAND' if mb_b[i]>mb_h[i] else 'HOP'
    print(f"  T={Tv:4.0f}K: Band={mb_b[i]:.3f} Hop={mb_h[i]:.3f} B/H={mb_b[i]/mb_h[i]:.2f} -> {d}")

save_checkpoint({**{f'mu_tot_{d}':r_mb[f'mu_total_{d}'] for d in ['a','b','c_prime']},
    **{f'mu_band_{d}':r_mb[f'mu_band_{d}'] for d in ['a','b','c_prime']},
    **{f'mu_hop_{d}':r_mb[f'mu_hopping_{d}'] for d in ['a','b','c_prime']},
    'T_range':Tr}, 'data/results/stage6_meanband.npz')
print("\nSaved.")
