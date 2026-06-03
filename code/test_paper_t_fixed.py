#!/usr/bin/env python3
"""Test paper t values with FIXED units."""
import sys; sys.path.insert(0,'.')
import numpy as np
from naph_transport.stage5_mobility_assemble import compute_direction_projections, compute_A0_alpha
from naph_transport.stage6_mobility_integrate import compute_mobility_vs_T
from naph_transport.io_utils import load_checkpoint
from naph_transport.constants import CM1_TO_RADS, REF_T_VALUES_EV

t_paper = np.array([REF_T_VALUES_EV[k] for k in ['a','b','c','ac','ab','abc']])
s5 = load_checkpoint('data/results/stage5_mobility_params.npz')
omega = s5['omega_cm1'] * CM1_TO_RADS; G = s5['G_total']
R_proj = compute_direction_projections()
A0_p = compute_A0_alpha(t_paper, R_proj)

T_range = np.arange(10,310,10)
print(f'{"T":>5s} {"mu_a":>8s} {"mu_b":>8s} {"mu_c":>8s}  {"band_c":>8s} {"hop_c":>8s}')
for d in ['a','b','c_prime']:
    r = compute_mobility_vs_T(T_range, A0_p[d], s5[f'B_{d}'], G, omega, verbose=False)
    locals()[f'r_{d}'] = r

for i,Tv in enumerate(T_range):
    if Tv in [10,50,100,200,300]:
        print(f'{Tv:5.0f} {r_a["mu_total"][i]:8.2f} {r_b["mu_total"][i]:8.2f} '
              f'{r_c_prime["mu_total"][i]:8.3f}  {r_c_prime["mu_band"][i]:8.3f} '
              f'{r_c_prime["mu_hopping"][i]:8.3f}')

mc = r_c_prime['mu_total'][-1]
print(f'\nmu_c(300K) paper t = {mc:.3f} cm2/Vs')
print(f'Band={r_c_prime["mu_band"][-1]:.3f}, Hop={r_c_prime["mu_hopping"][-1]:.3f}')
mechanism = 'Hopping-dominated' if r_c_prime['mu_hopping'][-1] > r_c_prime['mu_band'][-1] else 'Band-dominated'
print(f'Mechanism: {mechanism}')
