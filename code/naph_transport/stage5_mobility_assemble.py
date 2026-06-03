"""Stage 5: Mobility parameter assembly — A_alpha, B_alpha_q, G_lambda.

Reads Stage 1-4 checkpoints and computes band/hopping prefactors
for the Holstein-Peierls mobility formula:

  mu_alpha(T) = (e0/(kB*T)) * [A0_alpha * I0(T) + sum_lambda B_alpha_lambda * I_lambda(T)]

BUG FIXES from v1 audit:
  F1: B_alpha prefactor uses e0/(4*kB) — hbar cancels completely
  M11: Bose occupation clamped at exp arg < 700
"""
import numpy as np
from typing import Dict

from .constants import (E0_SI, KB_SI, HBAR_SI, EV_TO_J, A_TO_M, M2_TO_CM2,
                         CM1_TO_RADS, DIMER_PAIRS_SAME, DIMER_PAIRS_CROSS,
                         dE_dQ_to_dimensionless_g)
from .crystal_utils import CrystalParams
from .io_utils import load_checkpoint, save_checkpoint


def compute_direction_projections() -> Dict[str, np.ndarray]:
    """Compute R_alpha,n = |R_n · alpha_hat| for all 6 pair vectors.

    Returns:
        dict: {'a': (6,) array, 'b': (6,) array, 'c_prime': (6,) array}
    """
    cp = CrystalParams()
    lat = cp.lattice_vectors()
    a_vec = lat[0]; b_vec = lat[1]
    a_hat = a_vec / np.linalg.norm(a_vec)
    b_hat = b_vec / np.linalg.norm(b_vec)
    ab_normal = np.cross(a_hat, b_hat)
    c_prime_hat = ab_normal / np.linalg.norm(ab_normal)

    R_vectors = np.zeros((6, 3))
    labels = ['a', 'b', 'c', 'ac', 'ab', 'abc']
    for i, key in enumerate(labels[:3]):
        R_vectors[i] = cp.R_frac_to_cart(DIMER_PAIRS_SAME[key])
    for i, key in enumerate(labels[3:]):
        R_vectors[i + 3] = cp.R_frac_to_cart(DIMER_PAIRS_CROSS[key])

    projections = {}
    for name, hat in [('a', a_hat), ('b', b_hat), ('c_prime', c_prime_hat)]:
        projections[name] = np.abs(np.dot(R_vectors, hat))
    return projections


def compute_A0_alpha(t_values: np.ndarray, R_proj: Dict[str, np.ndarray]
                     ) -> Dict[str, float]:
    """Compute band prefactor A0_alpha for each direction.

    A0_alpha = (1/(2*hbar^2)) * sum_n R_alpha,n^2 * t_n^2

    Args:
        t_values: (6,) transfer integrals in eV [t_a, t_b, t_c, t_ac, t_ab, t_abc]
        R_proj: direction projections from compute_direction_projections()
    """
    A0 = {}
    t_si = t_values * EV_TO_J
    for direction in ['a', 'b', 'c_prime']:
        R_m = R_proj[direction] * A_TO_M
        A0[direction] = np.sum(R_m**2 * t_si**2) / (2.0 * HBAR_SI**2)
    return A0


def compute_G_total(g_local: np.ndarray, g_nonlocal: np.ndarray) -> np.ndarray:
    """G_lambda = |g_local|^2 + sum_n |g_n_nonlocal|^2."""
    return g_local**2 + np.sum(g_nonlocal**2, axis=1)


def compute_B_alpha_lambda(g_nonlocal: np.ndarray, omega_cm1: np.ndarray,
                            R_proj: Dict[str, np.ndarray]
                            ) -> Dict[str, np.ndarray]:
    """Compute mode-resolved hopping prefactor B_alpha_lambda.

    B_alpha_lambda = (omega_lambda^2 / 4) * sum_n R_alpha,n^2 * |g_n|^2
    F1 FIX: hbar cancels completely.
    """
    omega_rads = omega_cm1 * CM1_TO_RADS
    n_modes = len(omega_rads)
    B = {}
    for direction in ['a', 'b', 'c_prime']:
        R_m = R_proj[direction] * A_TO_M
        B_dir = np.zeros(n_modes)
        for lam in range(n_modes):
            B_dir[lam] = np.sum(R_m**2 * g_nonlocal[lam, :]**2)
        B_dir *= omega_rads**2 / 4.0
        B[direction] = B_dir
    return B


def run_stage5(
    stage1_path: str = 'data/results/stage1_transfer_integrals.npz',
    stage3_path: str = 'data/results/stage3_eph_coupling.npz',
    stage2_path: str = 'data/results/stage2_phonon.npz',
    output_path: str = 'data/results/stage5_mobility_params.npz',
    verbose: bool = True
) -> Dict:
    """Run Stage 5: assemble A0_alpha and B_alpha_lambda."""
    s1 = load_checkpoint(stage1_path)
    s2 = load_checkpoint(stage2_path)
    s3 = load_checkpoint(stage3_path)

    t_values = np.array([
        float(s1['t_a']), float(s1['t_b']), float(s1['t_c']),
        float(s1['t_ac']), float(s1['t_ab']), float(s1['t_abc'])
    ])
    omega_cm1 = s2['selected_frequencies']
    d_eps_dQ = s3['d_epsilon0_dQ']
    dt_dQ = s3['dt_dQ']

    n_modes = len(omega_cm1)
    g_local = np.array([dE_dQ_to_dimensionless_g(d_eps_dQ[i], abs(omega_cm1[i]))
                         for i in range(n_modes)])
    g_nonlocal = np.array([[dE_dQ_to_dimensionless_g(dt_dQ[i, j], abs(omega_cm1[i]))
                              for j in range(6)] for i in range(n_modes)])

    R_proj = compute_direction_projections()
    A0 = compute_A0_alpha(t_values, R_proj)
    G_total = compute_G_total(g_local, g_nonlocal)
    B = compute_B_alpha_lambda(g_nonlocal, omega_cm1, R_proj)

    # Keep A0 in SI units (m2/s2) — M2_TO_CM2 applied once in stage6
    if verbose:
        print("Stage 5: Mobility Parameter Assembly")
        print(f"  t (meV): {[f'{t*1000:.1f}' for t in t_values]}")
        print(f"  Modes: {n_modes}, ω range: [{omega_cm1.min():.0f}, {omega_cm1.max():.0f}] cm⁻¹")
        print(f"  G_total: [{G_total.min():.4f}, {G_total.max():.4f}]")
        for d in ['a', 'b', 'c_prime']:
            print(f"  A0_{d}: {A0[d]*M2_TO_CM2:.2e} cm2/s2 (SI→cm2), max|B_{d}|: {np.max(np.abs(B[d])):.2e}")

    save_checkpoint({
        't_values': t_values, 'omega_cm1': omega_cm1,
        'G_total': G_total, 'g_local': g_local, 'g_nonlocal': g_nonlocal,
        'A0_a': A0['a'], 'A0_b': A0['b'], 'A0_c_prime': A0['c_prime'],
        'B_a': B['a'], 'B_b': B['b'], 'B_c_prime': B['c_prime'],
    }, output_path)

    return {'A0': A0, 'B': B, 'G_total': G_total, 'R_proj': R_proj}
