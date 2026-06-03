"""Validation checkpoints for the naphthalene transport pipeline.

CP1: TB fitting quality (condition number, RMSE, parameter comparison)
CP2: Phonon spectrum quality (acoustic modes, no imaginary frequencies)
CP3: E-ph coupling linearity (parameter antisymmetry, force antisymmetry)
CP4: Assembled parameter magnitudes
"""
import numpy as np
from typing import Dict


def cp1_check_tb_fit(condition_number: float, rmse_eV: float,
                     t_fitted: Dict[str, float],
                     t_ref: Dict[str, float],
                     verbose: bool = True) -> Dict[str, bool]:
    """Checkpoint 1: Validate TB fitting results.

    Args:
        condition_number: design matrix condition number κ
        rmse_eV: root-mean-square error (eV)
        t_fitted: fitted transfer integrals {'a': t_a, 'b': t_b, ...}
        t_ref: reference transfer integrals from paper Table I
        verbose: print detailed results

    Returns:
        dict of check_name → passed (bool)
    """
    results = {}

    # Check 1: Condition number
    results['cond_ok'] = condition_number < 100.0
    if verbose:
        print(f"CP1.1 Condition number: κ={condition_number:.1f} "
              f"({'PASS' if results['cond_ok'] else 'FAIL — κ≥100'})")

    # Check 2: RMSE
    rmse_meV = rmse_eV * 1000.0
    results['rmse_ok'] = rmse_meV < 5.0
    if verbose:
        print(f"CP1.2 Fitting RMSE: {rmse_meV:.2f} meV "
              f"({'PASS' if results['rmse_ok'] else 'FAIL — RMSE≥5 meV'})")

    # Check 3: Parameter sign consistency
    signs_match = True
    for key in ['a', 'b', 'c', 'ac', 'ab', 'abc']:
        if key in t_fitted and key in t_ref:
            if np.sign(t_fitted[key]) != np.sign(t_ref[key]):
                signs_match = False
                if verbose:
                    print(f"  Sign mismatch for t_{key}: fitted={t_fitted[key]:.4f}, "
                          f"ref={t_ref[key]:.4f}")
    results['signs_ok'] = signs_match
    if verbose:
        print(f"CP1.3 Parameter signs: "
              f"{'PASS' if results['signs_ok'] else 'FAIL — sign mismatch'}")

    # Check 4: Parameter magnitude
    mag_ok = True
    for key in ['a', 'b', 'c', 'ac', 'ab', 'abc']:
        if key in t_fitted and key in t_ref:
            t_fit_abs = abs(t_fitted[key])
            t_ref_abs = abs(t_ref[key])
            if t_ref_abs > 1e-6:  # avoid division by near-zero reference
                ratio = t_fit_abs / t_ref_abs
                if ratio < 0.5 or ratio > 2.0:
                    mag_ok = False
                    if verbose:
                        print(f"  Magnitude mismatch for t_{key}: "
                              f"fitted={t_fitted[key]:.4f}, ref={t_ref[key]:.4f} "
                              f"(ratio={ratio:.2f})")
    results['magnitude_ok'] = mag_ok
    if verbose:
        print(f"CP1.4 Parameter magnitude: "
              f"{'PASS' if results['magnitude_ok'] else 'FAIL — >2× deviation'}")

    return results


def cp2_check_phonon(frequencies_cm1: np.ndarray,
                     verbose: bool = True) -> Dict[str, bool]:
    """Checkpoint 2: Validate phonon spectrum.

    Args:
        frequencies_cm1: (N_modes,) frequencies in cm⁻¹
        verbose: print detailed results

    Returns:
        dict of check_name → passed (bool)
    """
    results = {}

    # Check 1: Acoustic modes near zero (allow small DFPT truncation errors ±5 cm⁻¹)
    acoustic = np.sort(np.abs(frequencies_cm1[:3]))
    results['acoustic_ok'] = np.all(acoustic < 5.0)
    if verbose:
        print(f"CP2.1 Acoustic modes: {acoustic} cm⁻¹ "
              f"({'PASS' if results['acoustic_ok'] else 'FAIL — >5 cm⁻¹'})")

    # Check 2: No imaginary frequencies in optical modes (ω > +5 cm⁻¹)
    optical = frequencies_cm1[3:]
    imag_optical = np.sum(optical < -5.0)
    results['no_imag_optical'] = imag_optical == 0
    if verbose:
        print(f"CP2.2 Imaginary optical modes: {imag_optical} "
              f"({'PASS' if results['no_imag_optical'] else 'FAIL — unstable structure'})")

    # Check 3: C-H stretch frequencies in expected range
    max_freq = np.max(frequencies_cm1)
    results['ch_stretch_ok'] = 2800.0 <= max_freq <= 3500.0
    if verbose:
        print(f"CP2.3 Max frequency (C-H stretch): {max_freq:.1f} cm⁻¹ "
              f"({'PASS' if results['ch_stretch_ok'] else 'FAIL — out of 2800-3500 range'})")

    # Check 4: Low-frequency intermolecular modes exist
    intermolecular_count = np.sum((frequencies_cm1 > 30.0) & (frequencies_cm1 < 200.0))
    results['intermol_ok'] = intermolecular_count >= 5
    if verbose:
        print(f"CP2.4 Intermolecular modes (30-200 cm⁻¹): {intermolecular_count} "
              f"({'PASS' if results['intermol_ok'] else 'FAIL — <5 modes'})")

    return results


def cp3_check_eph_linearity(t_plus: Dict[str, float], t_zero: Dict[str, float],
                             t_minus: Dict[str, float], forces_plus: np.ndarray,
                             forces_minus: np.ndarray, verbose: bool = True
                             ) -> Dict[str, bool]:
    """Checkpoint 3: Validate e-ph coupling linearity.

    Checks antisymmetry of t parameters and force antisymmetry for a single mode.

    Args:
        t_plus: TB parameters at +ΔQ
        t_zero: TB parameters at equilibrium
        t_minus: TB parameters at -ΔQ
        forces_plus: (N_atoms, 3) forces at +ΔQ (eV/Å)
        forces_minus: (N_atoms, 3) forces at -ΔQ (eV/Å)
        verbose: print detailed results

    Returns:
        dict of check_name → passed (bool)
    """
    results = {}

    # Check 1: Parameter antisymmetry (t(+ΔQ) - t(0) ≈ t(0) - t(-ΔQ))
    max_delta_t = 0.0
    for key in t_plus:
        if key in t_zero and key in t_minus:
            forward = abs(t_plus[key] - t_zero[key])
            backward = abs(t_zero[key] - t_minus[key])
            if forward + backward > 1e-10:
                asym = abs(forward - backward) / (forward + backward)
                max_delta_t = max(max_delta_t, asym)

    results['t_symmetry_ok'] = max_delta_t < 0.2  # 20% asymmetry tolerance
    if verbose:
        print(f"CP3.1 TB parameter antisymmetry (max asymmetry): {max_delta_t:.4f} "
              f"({'PASS' if results['t_symmetry_ok'] else 'FAIL — >20% asymmetric'})")

    # Check 2: Force antisymmetry (F⁺ + F⁻ ≈ 0)
    force_sum = np.max(np.abs(forces_plus + forces_minus))
    results['force_antisym_ok'] = force_sum < 0.01  # eV/Å
    if verbose:
        print(f"CP3.2 Force antisymmetry max|F⁺+F⁻|: {force_sum:.6f} eV/Å "
              f"({'PASS' if results['force_antisym_ok'] else 'FAIL — >0.01 eV/Å'})")

    return results


def cp4_check_assembled(A_alpha: Dict[str, float],
                        G_lambda: np.ndarray,
                        omega_cm1: np.ndarray,
                        verbose: bool = True) -> Dict[str, bool]:
    """Checkpoint 4: Validate assembled parameters.

    Args:
        A_alpha: band prefactors {'a': A_a, 'b': A_b, 'c': A_c} in cm²·K/V·s²
        G_lambda: (N_modes,) dimensionless coupling constants (HOMO)
        omega_cm1: (N_modes,) mode frequencies in cm⁻¹
        verbose: print detailed results

    Returns:
        dict of check_name → passed (bool)
    """
    results = {}

    # Check 1: A_alpha magnitude
    A_values = list(A_alpha.values())
    A_ok = all(1e14 < abs(a) < 1e19 for a in A_values)
    results['A_mag_ok'] = A_ok
    if verbose:
        for key, val in A_alpha.items():
            print(f"CP4.1 A_{key} = {val:.2e} cm²·K/V·s²")
        print(f"CP4.1 A_alpha magnitude: {'PASS' if A_ok else 'FAIL'}")

    # Check 2: G_lambda magnitude for HOMO (should be < 1 for most modes)
    max_G = np.max(np.abs(G_lambda))
    results['G_mag_ok'] = max_G < 1.5
    if verbose:
        print(f"CP4.2 Max |G_λ| (HOMO): {max_G:.4f} "
              f"({'PASS' if results['G_mag_ok'] else 'FAIL — >1.5'})")

    # Check 3: Low-frequency modes dominate
    low_mask = omega_cm1 < 200.0
    if np.sum(low_mask) > 0:
        G_low_mean = np.mean(np.abs(G_lambda[low_mask]))
        G_high_mean = np.mean(np.abs(G_lambda[~low_mask])) if np.sum(~low_mask) > 0 else 0.0
        results['low_freq_dom_ok'] = G_low_mean > G_high_mean
        if verbose:
            print(f"CP4.3 Mean |G| low-freq (<200 cm⁻¹): {G_low_mean:.4f}, "
                  f"high-freq: {G_high_mean:.4f} "
                  f"({'PASS' if results['low_freq_dom_ok'] else 'FAIL'})")
    else:
        results['low_freq_dom_ok'] = True

    return results
