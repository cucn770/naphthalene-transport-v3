"""Stage 6: Mobility via time-domain numerical integration.

Holstein-Peierls model (J. Chem. Phys. 127, 044506, 2007, Eq. 11-15):

  mu_alpha(T) = (e0/(kB*T)) * [A0_alpha * I0(T) + sum_lambda B_alpha_lambda * I_lambda(T)]

BUG FIXES from v1: F2 adaptive cutoff, S8 Simpson+Nyquist, M10 convergence, M11 Bose clamp.
"""
import numpy as np
from scipy.integrate import simpson
from typing import Tuple, Dict

from .constants import (E0_SI, KB_SI, HBAR_SI, HBAR_EV, KB_EV, CM1_TO_RADS, M2_TO_CM2)
from .io_utils import load_checkpoint, save_checkpoint


def bose_occupation(omega_rads: np.ndarray, T: float) -> np.ndarray:
    """Bose-Einstein N_lambda(T). M11: clamp exp arg < 700."""
    exponent = omega_rads * HBAR_EV / (KB_EV * T)
    N = np.zeros_like(exponent)
    mask = exponent < 700.0
    N[mask] = 1.0 / (np.exp(exponent[mask]) - 1.0)
    return N


def determine_integration_params(G_total, N_T, omega_rads, Gamma_eV=1e-4):
    """F2: Adaptive cutoff from Gaussian + Debye-Waller decay."""
    Gamma_rad = Gamma_eV / HBAR_EV
    t_cut_gauss = 5.0 / Gamma_rad if Gamma_rad > 0 else np.inf
    dw_coeff = np.sum(G_total * (1.0 + 2.0 * N_T) * omega_rads**2)
    t_cut_dw = 5.0 / np.sqrt(dw_coeff) if dw_coeff > 0 else np.inf
    t_cut = min(t_cut_gauss, t_cut_dw, 1e-12)
    omega_max = np.max(omega_rads[omega_rads > 0]) if np.any(omega_rads > 0) else 1e14
    dt = min(0.5 * np.pi / omega_max, 2e-15)
    dt = max(dt, 1e-17)
    n_points = max(int(t_cut / dt) + 1, 100)
    return t_cut, dt, n_points


def compute_C_and_S(t_array, G_total, N_T, omega_rads, Gamma_rad):
    """C(t)=exp[-2ΣG(1+2N)(1-cos ωt)]·exp[-Γ²t²], S(t)=2ΣG sin ωt."""
    n_t = len(t_array); n_modes = len(G_total)
    sum_cos = np.zeros(n_t); sum_sin = np.zeros(n_t)
    for lam in range(n_modes):
        w_t = omega_rads[lam] * t_array
        pre = 2.0 * G_total[lam] * (1.0 + 2.0 * N_T[lam])
        sum_cos += pre * (1.0 - np.cos(w_t))
        sum_sin += 2.0 * G_total[lam] * np.sin(w_t)
    return (np.exp(-sum_cos) * np.exp(-Gamma_rad**2 * t_array**2),
            sum_sin)


def compute_mobility_at_T(T, A0_alpha, B_alpha, G_total, omega_rads,
                           Gamma_eV=1e-4, verbose=False):
    """μ = (e0/(kB*T)) * [A0*I0 + Σ B_λ*I_λ]. Returns (total, band, hopping)."""
    N_T = bose_occupation(omega_rads, T)
    t_cut, dt, n_points = determine_integration_params(G_total, N_T, omega_rads, Gamma_eV)
    t_array = np.linspace(0, t_cut, n_points)
    Gamma_rad = Gamma_eV / HBAR_EV
    C, S = compute_C_and_S(t_array, G_total, N_T, omega_rads, Gamma_rad)

    I0 = simpson(C * np.cos(S), x=t_array)

    hopping_sum = 0.0; n_modes = len(G_total)
    for lam in range(n_modes):
        if abs(omega_rads[lam]) < 1e-6: continue
        w_t = omega_rads[lam] * t_array
        I_lam = simpson(C * ((1+N_T[lam])*np.cos(S+w_t) + N_T[lam]*np.cos(S-w_t)), x=t_array)
        hopping_sum += B_alpha[lam] * I_lam

    prefactor = E0_SI / (KB_SI * T)
    mu_band = prefactor * A0_alpha * I0 * M2_TO_CM2
    mu_hopping = prefactor * hopping_sum * M2_TO_CM2
    mu_total = mu_band + mu_hopping

    if verbose and (abs(T-10)<1 or abs(T-300)<1 or abs(T%50)<1):
        print(f"  T={T:4.0f}K μ={mu_total:.4f} (band={mu_band:.4f} hop={mu_hopping:.4f}) cm²/Vs")
    return mu_total, mu_band, mu_hopping


def compute_mobility_vs_T(T_range, A0_alpha, B_alpha, G_total, omega_rads,
                           Gamma_eV=1e-4, verbose=True):
    """Compute μ(T) for all temperatures."""
    n_T = len(T_range)
    mu_tot, mu_band, mu_hop = np.zeros(n_T), np.zeros(n_T), np.zeros(n_T)
    for i, T in enumerate(T_range):
        mu_tot[i], mu_band[i], mu_hop[i] = compute_mobility_at_T(
            T, A0_alpha, B_alpha, G_total, omega_rads, Gamma_eV, verbose=verbose)
    return {'mu_total': mu_tot, 'mu_band': mu_band, 'mu_hopping': mu_hop}


def verify_convergence(T, A0, B, G, omega, Gamma_eV=1e-4):
    """M10: dt-halving convergence test."""
    mu_ref, _, _ = compute_mobility_at_T(T, A0, B, G, omega, Gamma_eV)
    N_T = bose_occupation(omega, T)
    t_cut, dt, n_pts = determine_integration_params(G, N_T, omega, Gamma_eV)
    t_fine = np.linspace(0, t_cut, n_pts * 2)
    Gamma_rad = Gamma_eV / HBAR_EV
    C, S = compute_C_and_S(t_fine, G, N_T, omega, Gamma_rad)
    I0f = simpson(C * np.cos(S), x=t_fine)
    hs = 0.0
    for lam in range(len(G)):
        if abs(omega[lam]) < 1e-6: continue
        w_t = omega[lam] * t_fine
        hs += B[lam] * simpson(C * ((1+N_T[lam])*np.cos(S+w_t)+N_T[lam]*np.cos(S-w_t)), x=t_fine)
    mu_fine = (E0_SI/(KB_SI*T)) * (A0*I0f + hs) * M2_TO_CM2
    delta = abs(mu_fine - mu_ref) / (abs(mu_ref) + 1e-30)
    ok = delta < 0.05
    print(f"  Convergence T={T}K: μ={mu_ref:.4f} μ_fine={mu_fine:.4f} Δ={delta*100:.1f}% "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def run_stage6(stage5_path='data/results/stage5_mobility_params.npz',
               output_path='data/results/stage6_mobility_results.npz',
               T_min=10.0, T_max=300.0, T_step=10.0,
               Gamma_eV=1e-4, verbose=True):
    """Run Stage 6: μ(T) for a, b, c' directions."""
    s5 = load_checkpoint(stage5_path)
    omega_rads = s5['omega_cm1'] * CM1_TO_RADS
    G_total = s5['G_total']
    T_range = np.arange(T_min, T_max + T_step, T_step)
    results = {'T_range': T_range}

    if verbose:
        print(f"Stage 6: Mobility ({len(G_total)} modes, Γ={Gamma_eV*1000:.1f} meV)")

    for d in ['a', 'b', 'c_prime']:
        A0 = float(s5[f'A0_{d}'])
        B = s5[f'B_{d}']
        if verbose: print(f"\n  Direction {d}: A0={A0:.2e}")
        dr = compute_mobility_vs_T(T_range, A0, B, G_total, omega_rads, Gamma_eV, verbose)
        for k, v in dr.items(): results[f'{k}_{d}'] = v
        if verbose: verify_convergence(T_max, A0, B, G_total, omega_rads, Gamma_eV)

    save_checkpoint(results, output_path)
    if verbose:
        print(f"\nStage 6 done. μ(300K):")
        for d in ['a','b','c_prime']: print(f"  {d}: {results[f'mu_total_{d}'][-1]:.4f} cm²/Vs")
    return results
