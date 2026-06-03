"""Stage 3: Electron-phonon coupling from finite-difference TB fitting.

For each phonon mode λ:
  1. Generate ±ΔQ displaced structures (via Phonopy)
  2. Run VASP single-point band calculation on displaced structures
  3. Re-fit 7-parameter TB model to extract ε₀(±ΔQ), t_mn(±ΔQ)
  4. Central finite difference: ∂ε₀/∂Q, ∂t_mn/∂Q

Key: ∂ε₀/∂Q gives the pure Holstein (local) coupling.
"""
import numpy as np
import os
from typing import Dict, Tuple, List

from .io_utils import (parse_vasp_eigenval, parse_vasp_outcar_forces,
                        save_checkpoint, load_checkpoint)
from .stage1_transfer import fit_tb_parameters
from .constants import DIMER_PAIRS_SAME, DIMER_PAIRS_CROSS, REF_T_VALUES_EV
from .validator import cp3_check_eph_linearity


def enforce_sign_continuity(t_fitted: Dict[str, float],
                             t_equilibrium: Dict[str, float],
                             mode_idx: int = 0,
                             verbose: bool = True
                             ) -> Tuple[Dict[str, float], list]:
    """Enforce sign continuity for cross-sublattice TB parameters.

    For small displacements (Delta_Q = 0.01 sqrt(amu)*Ang), transfer
    integrals should not flip sign. If sign changes but magnitude change
    is small (< 50%), it's a fitting artifact -- correct the sign.
    If magnitude change is large, warn for manual inspection.

    Args:
        t_fitted: TB parameters at +/- Delta_Q
        t_equilibrium: TB parameters at equilibrium (Stage 1 v3 result)
        mode_idx: mode index for logging
        verbose: print corrections

    Returns:
        t_fixed: sign-corrected parameters
        warnings: list of warning strings for large sign flips
    """
    t_fixed = dict(t_fitted)
    warnings = []
    # Only check cross-sublattice parameters (sign-sensitive in eigenvalue fitting)
    cross_keys = ['ac', 'ab', 'abc']
    for key in cross_keys:
        if key not in t_fitted or key not in t_equilibrium:
            continue
        t_eq = t_equilibrium[key]
        t_fit = t_fitted[key]

        if np.sign(t_fit) != np.sign(t_eq):
            if abs(t_fit) < 1e-8:
                # Near-zero: sign is numerically meaningless, keep as-is
                continue
            rel_change = abs(t_fit - t_eq) / max(abs(t_eq), 1e-6)
            if rel_change < 0.5:
                # Small relative change, sign flip is fitting artifact -> correct
                t_fixed[key] = np.sign(t_eq) * abs(t_fit)
                if verbose:
                    print(f"  [Mode {mode_idx}] {key}: sign corrected "
                          f"({t_fit*1000:.2f} -> {t_fixed[key]*1000:.2f} meV, "
                          f"delta_rel={rel_change*100:.0f}%)")
            else:
                # Large change -- could be real physics or fitting failure
                warnings.append(
                    f"Mode {mode_idx}, {key}: sign flip with large |Delta|t| "
                    f"({rel_change*100:.0f}%), t_eq={t_eq*1000:.1f}, "
                    f"t_fit={t_fit*1000:.1f} meV -- CHECK MANUALLY")
    return t_fixed, warnings


def process_one_displacement(eigenval_path: str, outcar_path: str,
                              tau_AB: np.ndarray,
                              t_init: Dict[str, float],
                              mode_idx: int = 0,
                              verbose: bool = False
                              ) -> Tuple[Dict[str, float], np.ndarray, list]:
    """Process one displaced structure: TB fit + force extraction + sign check.

    E4 FIX: Cross-sublattice transfer integrals are checked for sign
    continuity against the equilibrium result. Sign flips with small
    magnitude changes are corrected; large changes raise warnings.

    Args:
        eigenval_path: path to EIGENVAL
        outcar_path: path to OUTCAR
        tau_AB: sublattice displacement (Å)
        t_init: initial TB parameter guess (equilibrium values from Stage 1)
        mode_idx: mode index for logging
        verbose: print fitting progress

    Returns:
        t_fitted: sign-continuity-corrected TB parameters
        forces: (N_atoms, 3) forces in eV/Å
        sign_warnings: list of warning strings (empty if no issues)
    """
    k_frac, bands, n_elec = parse_vasp_eigenval(eigenval_path)
    homo_idx = n_elec // 2 - 1

    E_homo = bands[:, homo_idx]
    E_homo_minus_1 = bands[:, homo_idx - 1]

    t_fitted, _, _, _ = fit_tb_parameters(
        k_frac, E_homo, E_homo_minus_1, tau_AB, t_init=t_init, verbose=verbose)

    # E4 FIX: enforce sign continuity for cross-sublattice parameters
    t_fitted, sign_warnings = enforce_sign_continuity(
        t_fitted, t_init, mode_idx=mode_idx, verbose=verbose)

    forces = parse_vasp_outcar_forces(outcar_path)

    return t_fitted, forces, sign_warnings


def compute_dQ_derivatives(t_plus: Dict[str, float],
                            t_minus: Dict[str, float],
                            delta_Q: float = 0.01
                            ) -> Dict[str, float]:
    """Central finite difference: ∂P/∂Q = [P(+ΔQ) - P(-ΔQ)] / (2ΔQ).

    Args:
        t_plus: TB parameters at +ΔQ
        t_minus: TB parameters at -ΔQ
        delta_Q: displacement amplitude in √amu·Å

    Returns:
        dP_dQ: derivatives for each TB parameter
    """
    dP_dQ = {}
    for key in t_plus:
        if key in t_minus:
            dP_dQ[f'd{key}_dQ'] = (t_plus[key] - t_minus[key]) / (2.0 * delta_Q)
    return dP_dQ


def run_stage3(eph_dirs: List[str], t_equilibrium: Dict[str, float],
               tau_AB: np.ndarray,
               delta_Q: float = 0.01, verbose: bool = True) -> Dict:
    """Run Stage 3: process all e-ph displacement results.

    Args:
        eph_dirs: list of paths to displacement directories (+ and - alternating)
        t_equilibrium: TB parameters at equilibrium (from Stage 1)
        tau_AB: sublattice displacement (Å)
        delta_Q: displacement amplitude in √amu·Å
        verbose: print progress

    Returns:
        dict with ∂ε₀/∂Q and ∂t_mn/∂Q arrays
    """
    n_modes = len(eph_dirs) // 2
    pair_names = ['a', 'b', 'c', 'ac', 'ab', 'abc']

    d_epsilon0_dQ = np.zeros(n_modes)
    dt_dQ = np.zeros((n_modes, 6))
    forces_all = np.zeros((n_modes, 2, 36, 3))
    all_warnings = []

    for i in range(n_modes):
        dir_plus = eph_dirs[2 * i]
        dir_minus = eph_dirs[2 * i + 1]

        # Process +ΔQ (E4: sign continuity enforced)
        t_plus, f_plus, w_plus = process_one_displacement(
            os.path.join(dir_plus, 'EIGENVAL'),
            os.path.join(dir_plus, 'OUTCAR'),
            tau_AB, t_equilibrium, mode_idx=i)

        # Process -ΔQ (E4: sign continuity enforced)
        t_minus, f_minus, w_minus = process_one_displacement(
            os.path.join(dir_minus, 'EIGENVAL'),
            os.path.join(dir_minus, 'OUTCAR'),
            tau_AB, t_equilibrium, mode_idx=i)

        all_warnings.extend(w_plus)
        all_warnings.extend(w_minus)

        # Compute derivatives
        derivs = compute_dQ_derivatives(t_plus, t_minus, delta_Q)

        d_epsilon0_dQ[i] = derivs.get('depsilon_0_dQ', 0.0)
        for j, name in enumerate(pair_names):
            dt_dQ[i, j] = derivs.get(f'd{name}_dQ', 0.0)

        forces_all[i, 0] = f_plus
        forces_all[i, 1] = f_minus

        # Validate linearity for first few modes
        if verbose and i < 5:
            cp3_check_eph_linearity(t_plus, t_equilibrium, t_minus,
                                     f_plus, f_minus, verbose=True)

    if verbose:
        print(f"\nStage 3: E-ph coupling complete ({n_modes} modes)")
        if all_warnings:
            print(f"  WARNING: {len(all_warnings)} sign continuity issues:")
            for w in all_warnings:
                print(f"    {w}")
        print(f"  ∂ε₀/∂Q range: {np.min(d_epsilon0_dQ):.4f} to "
              f"{np.max(d_epsilon0_dQ):.4f} eV/√amu·Å")
        for j, name in enumerate(pair_names):
            print(f"  ∂t_{name}/∂Q range: {np.min(dt_dQ[:, j]):.4f} to "
                  f"{np.max(dt_dQ[:, j]):.4f} eV/√amu·Å")

    return {
        'd_epsilon0_dQ': d_epsilon0_dQ,
        'dt_dQ': dt_dQ,
        'forces': forces_all,
    }
