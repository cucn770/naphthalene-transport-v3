"""Stage 4: Assemble Holstein-Peierls model parameters.

Computes:
  - A_α: band transport prefactors (from t_mn)
  - G_λ: dimensionless e-ph coupling constants
  - B_αq: mode-resolved coupling prefactors

Does NOT compute mobility μ(T).
"""
import numpy as np
from typing import Dict

from .constants import (E0_SI, KB_SI, HBAR_SI, EV_TO_J, A_TO_M,
                         M2_TO_CM2, dE_dQ_to_dimensionless_g,
                         DIMER_PAIRS_SAME, DIMER_PAIRS_CROSS)
from .crystal_utils import CrystalParams
from .io_utils import save_checkpoint, load_checkpoint


def compute_G_dimensionless(d_epsilon0_dQ: np.ndarray,
                             d_t_dQ: np.ndarray,
                             frequencies_cm1: np.ndarray,
                             verbose: bool = True) -> Dict:
    """Compute dimensionless e-ph coupling constants G_λ.

    For each mode λ:
      g_λ^local = (∂ε₀/∂Q_λ) · ℓ_λ / (ħω_λ)          [Holstein]
      g_λ^nonlocal = (∂t_mn/∂Q_λ) · ℓ_λ / (ħω_λ)      [Peierls]

    where ℓ_λ = √(ħ/(2ω_λ)) is zero-point amplitude.

    Args:
        d_epsilon0_dQ: (N_modes,) ∂ε₀/∂Q in eV/√amu·Å
        d_t_dQ: (N_modes, 6) ∂t_mn/∂Q in eV/√amu·Å
        frequencies_cm1: (N_modes,) mode frequencies in cm⁻¹
        verbose: print results

    Returns:
        dict with G_local, G_nonlocal, G_total per mode
    """
    n_modes = len(frequencies_cm1)
    pair_names = ['a', 'b', 'c', 'ac', 'ab', 'abc']

    G_local = np.zeros(n_modes)
    G_nonlocal = np.zeros((n_modes, 6))

    for i in range(n_modes):
        omega_cm1 = abs(frequencies_cm1[i])
        if omega_cm1 < 1e-6:
            continue

        G_local[i] = dE_dQ_to_dimensionless_g(d_epsilon0_dQ[i], omega_cm1)

        for j in range(6):
            G_nonlocal[i, j] = dE_dQ_to_dimensionless_g(d_t_dQ[i, j], omega_cm1)

    G_total = np.sqrt(G_local**2 + np.sum(G_nonlocal**2, axis=1))

    if verbose:
        print(f"Stage 4: Parameter assembly")
        print(f"  G_local range: {np.min(np.abs(G_local)):.4f} to "
              f"{np.max(np.abs(G_local)):.4f}")
        print(f"  G_total range: {np.min(G_total):.4f} to {np.max(G_total):.4f}")
        top_idx = np.argsort(G_total)[::-1][:5]
        print(f"  Top 5 modes by G_total:")
        for idx in top_idx:
            print(f"    mode {idx+1}: ω={frequencies_cm1[idx]:.1f} cm⁻¹, "
                  f"G_local={G_local[idx]:.4f}, G_total={G_total[idx]:.4f}")

    return {
        'G_local': G_local,
        'G_nonlocal': G_nonlocal,
        'G_total': G_total,
        'pair_names': pair_names,
    }


def compute_A_alpha(t_equilibrium: Dict[str, float], verbose: bool = True
                    ) -> Dict[str, float]:
    """Compute band transport prefactors A_α.

    A_α ∝ t² · R² (dimensional prefactor).

    Args:
        t_equilibrium: fitted transfer integrals
        verbose: print results

    Returns:
        A_alpha in cm²·K/V·s²
    """
    cp = CrystalParams()
    A_prefactor_si = E0_SI / (2.0 * KB_SI * HBAR_SI**2)

    A_alpha = {}
    for direction, R_frac in {**DIMER_PAIRS_SAME, **DIMER_PAIRS_CROSS}.items():
        R_cart = cp.R_frac_to_cart(R_frac)
        R2 = np.dot(R_cart, R_cart)
        t = abs(t_equilibrium.get(direction, 0.0))
        t_si = t * EV_TO_J

        A_si = A_prefactor_si * R2 * (A_TO_M**2) * t_si**2
        A_alpha[direction] = A_si * M2_TO_CM2

    if verbose:
        for key, val in A_alpha.items():
            print(f"  A_{key} = {val:.2e} cm²·K/V·s²")

    return A_alpha


def run_stage4(stage1_path: str = 'data/results/stage1_transfer_integrals.npz',
               stage3_path: str = 'data/results/stage3_eph_coupling.npz',
               stage2_path: str = 'data/results/stage2_phonon.npz',
               output_path: str = 'data/results/stage4_assembled.npz',
               verbose: bool = True) -> Dict:
    """Run Stage 4: assemble all parameters.

    Args:
        stage1_path: Stage 1 checkpoint
        stage3_path: Stage 3 checkpoint
        stage2_path: Stage 2 checkpoint
        output_path: output checkpoint path
        verbose: print progress

    Returns:
        Assembled parameter dict
    """
    s1 = load_checkpoint(stage1_path)
    s2 = load_checkpoint(stage2_path)
    s3 = load_checkpoint(stage3_path)

    t_equilibrium = {
        'epsilon_0': float(s1['epsilon_0']),
        'a': float(s1['t_a']), 'b': float(s1['t_b']), 'c': float(s1['t_c']),
        'ac': float(s1['t_ac']), 'ab': float(s1['t_ab']), 'abc': float(s1['t_abc']),
    }

    frequencies = s2['selected_frequencies']
    d_epsilon0_dQ = s3['d_epsilon0_dQ']
    d_t_dQ = s3['dt_dQ']

    G = compute_G_dimensionless(d_epsilon0_dQ, d_t_dQ, frequencies, verbose=verbose)
    A_alpha = compute_A_alpha(t_equilibrium, verbose=verbose)

    save_data = {
        't_equilibrium_epsilon_0': t_equilibrium['epsilon_0'],
        **{f't_{k}': t_equilibrium[k] for k in ['a', 'b', 'c', 'ac', 'ab', 'abc']},
        'G_local': G['G_local'],
        'G_total': G['G_total'],
        'frequencies_cm1': frequencies,
        'd_epsilon0_dQ': d_epsilon0_dQ,
        'd_t_dQ': d_t_dQ,
        **{f'A_{k}': A_alpha.get(k, 0.0) for k in ['a', 'b', 'c', 'ac', 'ab', 'abc']},
    }
    save_checkpoint(save_data, output_path)

    if verbose:
        print(f"\nStage 4 complete. Results saved to: {output_path}")

    return save_data
