"""Stage 1: 2x2 Tight-Binding band fitting --> transfer integrals t_mn.

Paper-consistent method (J. Chem. Phys. 127, 044506 (2007)):
  - 2x2 TB Hamiltonian for 2 molecules/primitive cell
  - E±(k) = epsilon_0 + H_AA(k) +/- |H_AB(k)|
  - Simultaneous LSQ fit to HOMO and HOMO-1 bands
  - 6-direction mask: only a,b,c (same-sublattice) and ac,ab,abc (cross-sublattice)
"""
import numpy as np
from scipy.optimize import least_squares
from typing import Dict, Tuple, Optional

from .io_utils import parse_vasp_eigenval, parse_vasp_contcar, save_checkpoint
from .crystal_utils import (CrystalParams, compute_sublattice_displacement,
                              build_tb_H_AB)
from .constants import (DIMER_PAIRS_SAME, DIMER_PAIRS_CROSS, REF_T_VALUES_EV)
from .validator import cp1_check_tb_fit


def fit_tb_parameters(k_frac: np.ndarray, E_homo_dft: np.ndarray,
                      E_homo_minus_1_dft: np.ndarray,
                      tau_AB: np.ndarray,
                      t_init: Optional[Dict[str, float]] = None,
                      verbose: bool = True
                      ) -> Tuple[Dict[str, float], float, float, np.ndarray]:
    """Fit 7-parameter 2x2 TB model to DFT bands (paper-consistent method).

    Hamiltonian:
      H(k) = [[epsilon_0 + H_AA(k),     H_AB(k)      ],
              [H_AB*(k),                epsilon_0 + H_AA(k)]]

    H_AA(k) = 2*t_a*cos(k.R_a) + 2*t_b*cos(k.R_b) + 2*t_c*cos(k.R_c)
    H_AB(k) = t_ac*exp(ik.R_ac) + t_ab*exp(ik.R_ab) + t_abc*exp(ik.R_abc)

    Eigenvalues: E±(k) = epsilon_0 + H_AA(k) +/- |H_AB(k)|

    7 unknown parameters fitted to 72 data points (2 bands x 36 k-points).

    Args:
        k_frac: (N_k, 3) k-point coordinates (fractional reciprocal)
        E_homo_dft: (N_k,) HOMO band energies (eV)  --> E_+(k)
        E_homo_minus_1_dft: (N_k,) HOMO-1 energies (eV) --> E_-(k)
        tau_AB: (3,) sublattice displacement A-->B in Cartesian (Angstrom)
        t_init: initial guess (defaults to paper Table I)
        verbose: print fitting progress

    Returns:
        t_result: {'epsilon_0': e0, 'a': t_a, 'b': t_b, ...}
        rmse_eV: RMSE over 72 data points
        cond_num: Jacobian condition number
        E_fitted: (N_k, 2) fitted [E_plus, E_minus]
    """
    cp = CrystalParams()
    k_cart = cp.k_frac_to_cart(k_frac)

    # Build R vectors (Cartesian, Angstrom)
    R_same = np.array([cp.R_frac_to_cart(DIMER_PAIRS_SAME[k])
                       for k in ['a', 'b', 'c']])
    R_cross = np.array([cp.R_frac_to_cart(DIMER_PAIRS_CROSS[k])
                        for k in ['ac', 'ab', 'abc']])

    if verbose:
        print(f"2x2 TB model (paper-consistent method)")
        print(f"  tau_AB = ({tau_AB[0]:.3f}, {tau_AB[1]:.3f}, {tau_AB[2]:.3f}) A")
        for i, (label, R) in enumerate(zip(['a','b','c'], R_same)):
            print(f"  R_{label}: |R|={np.linalg.norm(R):.3f} A  [same-sublattice, H_AA]")
        for i, (label, R) in enumerate(zip(['ac','ab','abc'], R_cross)):
            print(f"  R_{label}: |R|={np.linalg.norm(R):.3f} A  [cross-sublattice, H_AB]")
        print(f"  E_HOMO range: [{E_homo_dft.min():.4f}, {E_homo_dft.max():.4f}] eV")
        print(f"  E_HOMO-1 range: [{E_homo_minus_1_dft.min():.4f}, {E_homo_minus_1_dft.max():.4f}] eV")
        split = (E_homo_dft - E_homo_minus_1_dft) / 2
        print(f"  Splitting: [{split.min()*1000:.1f}, {split.max()*1000:.1f}] meV")

    # Initial parameters
    if t_init is None:
        t_init = REF_T_VALUES_EV
    p0 = np.array([t_init['epsilon_0'],
                   t_init['a'], t_init['b'], t_init['c'],
                   t_init['ac'], t_init['ab'], t_init['abc']])

    def residuals(p):
        eps0, ta, tb, tc, tac, tab, tabc = p

        # Same-sublattice: real, cosine form
        H_AA = (2.0 * ta * np.cos(np.dot(k_cart, R_same[0])) +
                2.0 * tb * np.cos(np.dot(k_cart, R_same[1])) +
                2.0 * tc * np.cos(np.dot(k_cart, R_same[2])))

        # Cross-sublattice: P2_1/a symmetry via build_tb_H_AB
        # t_n * [exp(ik·(R+τ)) + exp(ik·(R-τ))] = 2·t_n·cos(k·τ)·exp(ik·R)
        t_cross_arr = np.array([tac, tab, tabc])
        H_AB = build_tb_H_AB(k_cart, R_cross, t_cross_arr, tau_AB)

        H_AB_abs = np.abs(H_AB)

        E_tb_plus = eps0 + H_AA + H_AB_abs
        E_tb_minus = eps0 + H_AA - H_AB_abs

        # Stack 72 residuals: [E_+(k0)-E_HOMO(k0), ..., E_-(k0)-E_HOMO-1(k0), ...]
        return np.concatenate([E_tb_plus - E_homo_dft,
                               E_tb_minus - E_homo_minus_1_dft])

    # Trust Region Reflective (well-posed problem with 6-direction mask)
    result = least_squares(residuals, p0, method='trf',
                           xtol=1e-14, ftol=1e-14, max_nfev=500,
                           verbose=0 if not verbose else 2)

    # Extract parameters
    labels = ['epsilon_0', 'a', 'b', 'c', 'ac', 'ab', 'abc']
    t_result = {lab: result.x[i] for i, lab in enumerate(labels)}

    # RMSE over all 72 data points
    rmse_eV = np.sqrt(np.mean(result.fun**2))

    # Condition number from final Jacobian
    J = result.jac
    _, S, _ = np.linalg.svd(J, full_matrices=False)
    cond_num = S.max() / S.min()

    # Compute fitted bands
    H_AA = (2.0 * t_result['a'] * np.cos(np.dot(k_cart, R_same[0])) +
            2.0 * t_result['b'] * np.cos(np.dot(k_cart, R_same[1])) +
            2.0 * t_result['c'] * np.cos(np.dot(k_cart, R_same[2])))
    t_cross_final = np.array([t_result['ac'], t_result['ab'], t_result['abc']])
    H_AB = build_tb_H_AB(k_cart, R_cross, t_cross_final, tau_AB)
    E_fitted_plus = t_result['epsilon_0'] + H_AA + np.abs(H_AB)
    E_fitted_minus = t_result['epsilon_0'] + H_AA - np.abs(H_AB)
    E_fitted = np.column_stack([E_fitted_plus, E_fitted_minus])

    if verbose:
        print(f"\n2x2 TB Fitting Results (72 data points, 7 parameters):")
        print(f"  epsilon_0 = {t_result['epsilon_0']*1000:.1f} meV  (paper: {REF_T_VALUES_EV['epsilon_0']*1000:.1f})")
        for key in ['a', 'b', 'c', 'ac', 'ab', 'abc']:
            sign_ok = "same" if np.sign(t_result[key]) == np.sign(REF_T_VALUES_EV[key]) else "OPP"
            print(f"  t_{key:>3s} = {t_result[key]*1000:7.2f} meV  "
                  f"(paper: {REF_T_VALUES_EV[key]*1000:7.2f})  sign: {sign_ok}")
        print(f"  RMSE = {rmse_eV*1000:.2f} meV  (over 72 data points)")
        print(f"  kappa = {cond_num:.1f}")

    return t_result, rmse_eV, cond_num, E_fitted


def run_stage1(eigenval_path: str = 'data/crystal/band/EIGENVAL',
               contcar_path: str = 'data/crystal/CONTCAR',
               output_path: str = 'data/results/stage1_transfer_integrals.npz',
               verbose: bool = True) -> Dict[str, float]:
    """Run Stage 1: 2x2 TB fitting (paper-consistent method).

    Args:
        eigenval_path: EIGENVAL path
        contcar_path: CONTCAR path
        output_path: output .npz path
        verbose: print progress

    Returns:
        t_result: fitted transfer integrals
    """
    k_frac, bands, n_elec = parse_vasp_eigenval(eigenval_path)
    homo_idx = n_elec // 2 - 1

    E_homo = bands[:, homo_idx]              # upper band  E_+
    E_homo_minus_1 = bands[:, homo_idx - 1]  # lower band  E_-

    lattice, pos_frac, _ = parse_vasp_contcar(contcar_path)
    tau_AB = compute_sublattice_displacement(pos_frac, lattice)

    t_result, rmse_eV, cond_num, E_fitted = fit_tb_parameters(
        k_frac, E_homo, E_homo_minus_1, tau_AB, verbose=verbose)

    # Validate
    cp1_check_tb_fit(cond_num, rmse_eV, t_result, REF_T_VALUES_EV, verbose=True)

    # Save
    save_data = {
        'k_frac': k_frac,
        'E_homo_dft': E_homo,
        'E_homo_minus_1_dft': E_homo_minus_1,
        'E_fitted_plus': E_fitted[:, 0],
        'E_fitted_minus': E_fitted[:, 1],
        'epsilon_0': t_result['epsilon_0'],
        't_a': t_result['a'], 't_b': t_result['b'], 't_c': t_result['c'],
        't_ac': t_result['ac'], 't_ab': t_result['ab'], 't_abc': t_result['abc'],
        'rmse_eV': rmse_eV, 'cond_num': cond_num, 'tau_AB': tau_AB,
    }
    save_checkpoint(save_data, output_path)

    if verbose:
        print(f"\nStage 1 complete. 2x2 paper-consistent method.")
        print(f"  RMSE={rmse_eV*1000:.2f} meV, kappa={cond_num:.1f}")
        print(f"  Bandwidth captured: {E_fitted[:,0].max()-E_fitted[:,0].min():.3f} eV "
              f"(DFT: {E_homo.max()-E_homo.min():.3f} eV)")
        print(f"Results: {output_path}")

    return t_result


# === E3 备选路径 (若build_tb_H_AB后RMSE未改善或恶化) ===
# P2_1/a 对称性要求对交叉子晶格dimer对考虑所有4个符号变体：
#
#   t_ab 连接 (±1/2, ±1/2, 0): 4个等效路径
#   t_abc 连接 (±1/2, ±1/2, 1): 4个等效路径
#   t_ac 连接 (1, 0, 1) 及其对称等效: 2个等效路径
#
# 展开形式示例:
#   H_AB = tac * (exp(ik·R_ac) + exp(-ik·R_ac))
#        + tab * (exp(ik·R_ab_pp) + exp(ik·R_ab_pm)
#               + exp(ik·R_ab_mp) + exp(ik·R_ab_mm))
#        + tabc * (exp(ik·R_abc_pp) + exp(ik·R_abc_pm)
#                + exp(ik·R_abc_mp) + exp(ik·R_abc_mm))
#
# 其中:
#   R_ab_pp = tau_AB + (0.5, 0.5, 0),  R_ab_pm = tau_AB + (0.5, -0.5, 0)
#   R_ab_mp = -tau_AB + (0.5, 0.5, 0), R_ab_mm = -tau_AB + (0.5, -0.5, 0)
#   (abc类似, 将最后分量改为1)
#
# 当前build_tb_H_AB使用 ±τ_AB 但固定R符号,
# 对应: 2·t_n·cos(k·τ)·exp(ik·R) = t_n·[exp(ik·(R+τ)) + exp(ik·(R-τ))]
