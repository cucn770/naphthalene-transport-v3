#!/usr/bin/env python3
"""Plot μ(T) following paper Fig 3/4 style. Compare v2 vs paper t vs paper reference."""
import sys; sys.path.insert(0,'.')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter

from naph_transport.io_utils import load_checkpoint
from naph_transport.stage5_mobility_assemble import (compute_direction_projections,
    compute_A0_alpha)
from naph_transport.stage6_mobility_integrate import compute_mobility_vs_T
from naph_transport.constants import CM1_TO_RADS, REF_T_VALUES_EV

# ---- Load v2 results ----
s6 = load_checkpoint('data/results/stage6_mobility_results.npz')
s5 = load_checkpoint('data/results/stage5_mobility_params.npz')
T_v2 = s6['T_range']

# ---- Compute paper-t mobility ----
t_paper = np.array([REF_T_VALUES_EV[k] for k in ['a','b','c','ac','ab','abc']])
R_proj = compute_direction_projections()
A0_p = compute_A0_alpha(t_paper, R_proj)
omega = s5['omega_cm1'] * CM1_TO_RADS
G = s5['G_total']

results_pt = {}
for d in ['a','b','c_prime']:
    r = compute_mobility_vs_T(T_v2, A0_p[d], s5[f'B_{d}'], G, omega, verbose=False)
    for k, v in r.items():
        results_pt[f'{k}_{d}'] = v

# ---- Paper reference data (approximate from Fig 3) ----
T_paper = np.array([10, 25, 50, 75, 100, 150, 200, 250, 300])
mu_paper_c = np.array([0.005, 0.008, 0.012, 0.022, 0.035, 0.055, 0.075, 0.09, 0.10])

# ============================================================
# Figure 1: μ(T) — log-log, paper Fig 3 style
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Left: v2 t values
ax = axes[0]
directions = [('a', 'blue'), ('b', 'red'), ('c_prime', 'green')]
labels = {'a': "a", 'b': "b", 'c_prime': "c'"}

for d, color in directions:
    mu_total = s6[f'mu_total_{d}']
    mu_band = s6[f'mu_band_{d}']
    mu_hop = s6[f'mu_hopping_{d}']
    ax.loglog(T_v2, mu_total, '-', color=color, linewidth=2, label=f"μ_{labels[d]} (v2 t)")
    ax.loglog(T_v2, mu_hop, '--', color=color, linewidth=1, alpha=0.5, label=f"hop_{labels[d]}")

ax.loglog(T_paper, mu_paper_c, 'ko-', markersize=6, label='Paper (c\')', zorder=10)
ax.set_xlabel('Temperature (K)', fontsize=12)
ax.set_ylabel('Mobility (cm²/Vs)', fontsize=12)
ax.set_title('v2 Transfer Integrals (2×2 model)', fontsize=13)
ax.legend(fontsize=8, loc='lower left')
ax.grid(True, alpha=0.3, which='both')
ax.set_xlim(8, 350)

# Right: paper t values
ax = axes[1]
for d, color in directions:
    mu_total = results_pt[f'mu_total_{d}']
    mu_band = results_pt[f'mu_band_{d}']
    mu_hop = results_pt[f'mu_hopping_{d}']
    ax.loglog(T_v2, mu_total, '-', color=color, linewidth=2, label=f"μ_{labels[d]} (paper t)")
    ax.loglog(T_v2, mu_hop, '--', color=color, linewidth=1, alpha=0.5, label=f"hop_{labels[d]}")

ax.loglog(T_paper, mu_paper_c, 'ko-', markersize=6, label='Paper ref (c\')', zorder=10)
ax.set_xlabel('Temperature (K)', fontsize=12)
ax.set_title('Paper Transfer Integrals (Table I)', fontsize=13)
ax.legend(fontsize=8, loc='lower left')
ax.grid(True, alpha=0.3, which='both')
ax.set_xlim(8, 350)

plt.suptitle('Naphthalene Hole Mobility: v2 vs Paper', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('results/figures/mobility_vs_T.png', dpi=150, bbox_inches='tight')
print("Saved: results/figures/mobility_vs_T.png")

# ============================================================
# Figure 2: Direction anisotropy — paper Fig 4 style
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: anisotropy ratio vs T (v2 t)
ax = axes[0]
mu_a_v2 = s6['mu_total_a']
mu_b_v2 = s6['mu_total_b']
mu_c_v2 = s6['mu_total_c_prime']
ax.plot(T_v2, mu_a_v2/mu_c_v2, 'b-', label='μ_a/μ_c\'')
ax.plot(T_v2, mu_b_v2/mu_c_v2, 'r-', label='μ_b/μ_c\'')
ax.plot(T_v2, mu_a_v2/mu_b_v2, 'g-', label='μ_a/μ_b')
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Anisotropy Ratio')
ax.set_title('v2 t values')
ax.legend()
ax.grid(True, alpha=0.3)

# Right: anisotropy ratio vs T (paper t)
ax = axes[1]
mu_a_pt = results_pt['mu_total_a']
mu_b_pt = results_pt['mu_total_b']
mu_c_pt = results_pt['mu_total_c_prime']
ax.plot(T_v2, mu_a_pt/mu_c_pt, 'b-', label='μ_a/μ_c\'')
ax.plot(T_v2, mu_b_pt/mu_c_pt, 'r-', label='μ_b/μ_c\'')
ax.plot(T_v2, mu_a_pt/mu_b_pt, 'g-', label='μ_a/μ_b')
ax.set_xlabel('Temperature (K)')
ax.set_title('Paper t values')
ax.legend()
ax.grid(True, alpha=0.3)

plt.suptitle('Mobility Anisotropy', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('results/figures/anisotropy.png', dpi=150, bbox_inches='tight')
print("Saved: results/figures/anisotropy.png")

# ============================================================
# Figure 3: Band vs Hopping decomposition (paper t values)
# ============================================================
fig, ax = plt.subplots(figsize=(10, 6))

mu_total_c = results_pt['mu_total_c_prime']
mu_band_c = results_pt['mu_band_c_prime']
mu_hop_c = results_pt['mu_hopping_c_prime']

ax.loglog(T_v2, mu_total_c, 'k-', linewidth=2, label='Total μ_c\'')
ax.loglog(T_v2, mu_band_c, 'b--', linewidth=1.5, label='Band (coherent)')
ax.loglog(T_v2, mu_hop_c, 'r--', linewidth=1.5, label='Hopping (incoherent)')
ax.loglog(T_paper, mu_paper_c, 'ko', markersize=8, label='Paper reference', zorder=10)

# Mark crossover region
crossover_T = T_v2[np.argmin(np.abs(mu_band_c - mu_hop_c))]
ax.axvline(x=crossover_T, color='gray', linestyle=':', alpha=0.5)
ax.text(crossover_T+5, ax.get_ylim()[0]*2, f'Band≈Hop\n~{crossover_T:.0f}K', fontsize=9)

ax.set_xlabel('Temperature (K)', fontsize=12)
ax.set_ylabel('Mobility (cm²/Vs)', fontsize=12)
ax.set_title('Naphthalene c\' Mobility: Band vs Hopping (Paper t values)', fontsize=13)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, which='both')
ax.set_xlim(8, 350)

plt.tight_layout()
plt.savefig('results/figures/band_vs_hopping.png', dpi=150, bbox_inches='tight')
print("Saved: results/figures/band_vs_hopping.png")

# ============================================================
# Figure 4: μ(T) comparison table (numerical values)
# ============================================================
print("\nNumerical results:")
print(f"{'T(K)':>6s} {'v2-μ_a':>10s} {'v2-μ_b':>10s} {'v2-μ_c':>10s}  {'pT-μ_a':>10s} {'pT-μ_b':>10s} {'pT-μ_c':>10s}")
print("-"*76)
for i, Tv in enumerate(T_v2):
    if Tv in [10, 25, 50, 75, 100, 150, 200, 250, 300]:
        print(f"{Tv:6.0f} {s6['mu_total_a'][i]:10.3f} {s6['mu_total_b'][i]:10.3f} {s6['mu_total_c_prime'][i]:10.3f}  "
              f"{results_pt['mu_total_a'][i]:10.3f} {results_pt['mu_total_b'][i]:10.3f} {results_pt['mu_total_c_prime'][i]:10.3f}")

print("\nAll figures saved to results/figures/")
