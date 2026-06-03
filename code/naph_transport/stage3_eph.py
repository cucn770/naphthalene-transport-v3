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


def process_one_displacement(eigenval_path: str, outcar_path: str,
                              tau_AB: np.ndarray,
                              t_init: Dict[str, float],
                              verbose: bool = False
                              ) -> Tuple[Dict[str, float], np.ndarray]:
    """Process one displaced structure: TB fit + force extraction.

    Args:
        eigenval_path: path to EIGENVAL
        outcar_path: path to OUTCAR
        tau_AB: sublattice displacement (Å)
        t_init: initial TB parameter guess
        verbose: print fitting progress

    Returns:
        t_fitted: fitted TB parameters
        forces: (N_atoms, 3) forces in eV/Å
    """
    k_frac, bands, n_elec = parse_vasp_eigenval(eigenval_path)
    homo_idx = n_elec // 2 - 1

    E_homo = bands[:, homo_idx]
    E_homo_minus_1 = bands[:, homo_idx - 1]

    t_fitted, _, _, _ = fit_tb_parameters(
        k_frac, E_homo, E_homo_minus_1, tau_AB, t_init=t_init, verbose=verbose)

    forces = parse_vasp_outcar_forces(outcar_path)

    return t_fitted, forces


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

    for i in range(n_modes):
        dir_plus = eph_dirs[2 * i]
        dir_minus = eph_dirs[2 * i + 1]

        # Process +ΔQ
        t_plus, f_plus = process_one_displacement(
            os.path.join(dir_plus, 'EIGENVAL'),
            os.path.join(dir_plus, 'OUTCAR'),
            tau_AB, t_equilibrium)

        # Process -ΔQ
        t_minus, f_minus = process_one_displacement(
            os.path.join(dir_minus, 'EIGENVAL'),
            os.path.join(dir_minus, 'OUTCAR'),
            tau_AB, t_equilibrium)

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
