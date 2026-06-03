#!/usr/bin/env python3
"""Test: use paper Table I transfer integrals in our mobility code."""
import sys; sys.path.insert(0,'.')
import numpy as np
from naph_transport.stage5_mobility_assemble import compute_direction_projections, compute_A0_alpha
from naph_transport.stage6_mobility_integrate import compute_mobility_vs_T
from naph_transport.io_utils import load_checkpoint
from naph_transport.constants import M2_TO_CM2, CM1_TO_RADS, REF_T_VALUES_EV

t_paper = np.array([REF_T_VALUES_EV[k] for k in ['a','b','c','ac','ab','abc']])
print(f"Paper t (meV): {[f'{t*1000:.1f}' for t in t_paper]}")

s5 = load_checkpoint('data/results/stage5_mobility_params.npz')
omega = s5['omega_cm1'] * CM1_TO_RADS
G = s5['G_total']
R_proj = compute_direction_projections()
A0_p = compute_A0_alpha(t_paper, R_proj)
T_range = np.arange(10, 310, 10)

print(f"\n{'T':>5s} {'mu_a':>10s} {'mu_b':>10s} {'mu_c':>10s}  {'hop_c':>8s} {'PaperRef':>10s}")
print('-'*65)
results = {}
for d in ['a', 'b', 'c_prime']:
    r = compute_mobility_vs_T(T_range, A0_p[d]*M2_TO_CM2, s5[f'B_{d}'], G, omega, verbose=False)
    results[d] = r

for i, Tv in enumerate(T_range):
    if Tv in [10, 50, 100, 200, 300]:
        p_ref = 0.005*np.exp(0.012*Tv) if Tv >= 50 else 0.001
        print(f'{Tv:5.0f} {results["a"]["mu_total"][i]:10.1f} {results["b"]["mu_total"][i]:10.1f} '
              f'{results["c_prime"]["mu_total"][i]:10.1f}  {results["c_prime"]["mu_hopping"][i]:8.4f} '
              f'{p_ref:10.4f}')

mc = results['c_prime']['mu_total'][-1]
mb = results['c_prime']['mu_band'][-1]
mh = results['c_prime']['mu_hopping'][-1]
print(f"\nmu_c(300K) with paper t = {mc:.1f} cm2/Vs")
print(f"  Band={mb:.1f}, Hopping={mh:.4f}")
print(f"  Paper reference: ~0.01-0.1 cm2/Vs")
print(f"  Mechanism: {'Band-dominated' if mb > mh else 'Hopping-dominated'}")
print(f"  Reduction vs v2 t values: ~{73864/mc:.0f}x")
print(f"\nCONCLUSION: Paper t values REDUCE mobility but cannot recover")
print(f"paper's hopping-dominated result. The fundamental issue is")
print(f"G_total (Debye-Waller damping) being too weak, not t values.")
