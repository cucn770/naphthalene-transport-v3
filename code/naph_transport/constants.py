"""Physical constants (CODATA 2018) and unit conversion utilities.

All internal calculations use SI units. Practical units (eV, cm⁻¹, Å, amu)
are converted to SI at input and back at output.

Key additions over v1:
  - Dimer pair classification (same-sublattice vs cross-sublattice)
  - Dimensionless coupling conversion factors
"""
import numpy as np

# ======== Fundamental Constants (SI) ========
HBAR_SI = 1.054571817e-34      # J·s
KB_SI = 1.380649e-23            # J/K
E0_SI = 1.602176634e-19         # C (elementary charge)
C_LIGHT = 2.99792458e8          # m/s
AMU_SI = 1.660539067e-27        # kg
HBAR_EV = 6.582119569e-16       # eV·s
KB_EV = 8.617333262e-5          # eV/K

# ======== Conversion Factors ========
EV_TO_J = 1.602176634e-19
J_TO_EV = 6.241509074e18
CM1_TO_EV = 1.239841984e-4          # cm⁻¹ → eV
EV_TO_CM1 = 8065.544                # eV → cm⁻¹
CM1_TO_RADS = 2.0 * np.pi * 2.99792458e10  # cm⁻¹ → rad/s
A_TO_M = 1.0e-10
A_TO_CM = 1.0e-8
M2_TO_CM2 = 1.0e4
AMU_SI_SQRT = np.sqrt(AMU_SI)       # √kg
AMU_HALF_A_TO_SI = AMU_SI_SQRT * A_TO_M  # √(amu)·Å → √(kg)·m

# ======== Naphthalene Crystal Structure ========
# P2_1/a space group, 2 molecules per primitive cell
# Molecules labeled sublattice A (molecule 1 in cell) and B (molecule 2 in cell)

# Dimer pair classification (from paper Table I):
#   Same-sublattice (AA/BB): translation by a, b, c → H_AA in TB model
#   Cross-sublattice (AB): translation by a+c, a/2±b/2, a/2±b/2+c → H_AB in TB model
#
# Translation vectors in direct (fractional) coordinates:
DIMER_PAIRS_SAME = {   # H_AA (same sublattice): real, symmetric → 2t·cos(k·R)
    'a': np.array([1.0, 0.0, 0.0]),
    'b': np.array([0.0, 1.0, 0.0]),
    'c': np.array([0.0, 0.0, 1.0]),
}

DIMER_PAIRS_CROSS = {  # H_AB (cross sublattice): complex phase factor
    'ac':  np.array([1.0, 0.0, 1.0]),
    # a/2+b/2 and a/2+b/2+c have 4 sign variants due to P2_1/a symmetry
    'ab':  np.array([0.5, 0.5, 0.0]),
    'abc': np.array([0.5, 0.5, 1.0]),
}

# Paper Table I reference values (eV) — for initial guess and validation
# From VASP PBE plane-wave periodic TB band-fitting
REF_T_VALUES_EV = {
    'epsilon_0': -0.835,
    'a':   -0.023,
    'b':   -0.042,
    'c':   -0.003,
    'ac':  -0.001,
    'ab':   0.022,
    'abc': -0.005,
}

# Paper Table II: intermolecular phonon frequencies (cm⁻¹)
REF_INTERMOLECULAR_CM1 = [58.8, 82.0, 108.5]

# Atom masses (amu)
MASS_C = 12.0107
MASS_H = 1.00794

# Sublattice definition:
# Molecule A = atoms 0-17 (indices from CONTCAR)
# Molecule B = atoms 18-35
MOL_A_ATOMS = list(range(0, 18))
MOL_B_ATOMS = list(range(18, 36))


def compute_zero_point_amplitude(omega_cm1):
    """Compute zero-point vibrational amplitude in √(amu)·Å.

    Args:
        omega_cm1: mode frequency in cm⁻¹

    Returns:
        ℓ_λ in √(amu)·Å (= √(ħ/(2Mω)) in mass-weighted coordinates)
    """
    omega_si = omega_cm1 * CM1_TO_RADS  # rad/s
    # ℓ in SI: √(ħ/(2ω)) in √(J·s²) = √(kg)·m
    ell_si = np.sqrt(HBAR_SI / (2.0 * omega_si))
    # Convert to √(amu)·Å
    return ell_si / AMU_HALF_A_TO_SI


def dE_dQ_to_dimensionless_g(dE_dQ_eV_per_amu_half_A, omega_cm1):
    """Convert energy derivative ∂E/∂Q to dimensionless coupling g.

    g = (∂E/∂Q) · ℓ / (ħω)

    Args:
        dE_dQ_eV_per_amu_half_A: ∂E/∂Q in eV/(√(amu)·Å)
        omega_cm1: mode frequency in cm⁻¹

    Returns:
        dimensionless coupling constant g
    """
    omega_si = omega_cm1 * CM1_TO_RADS  # rad/s
    hbar_omega_si = HBAR_SI * omega_si  # J
    hbar_omega_eV = omega_cm1 * CM1_TO_EV  # eV

    # ∂E/∂Q in SI: J/(√(kg)·m)
    dE_dQ_si = dE_dQ_eV_per_amu_half_A * EV_TO_J / AMU_HALF_A_TO_SI

    # Zero-point amplitude in √(kg)·m
    ell_si = np.sqrt(HBAR_SI / (2.0 * omega_si))

    # g = (∂E/∂Q) · ℓ / (ħω)
    g = dE_dQ_si * ell_si / hbar_omega_si
    return g


def cm1_to_rads(nu_cm1):
    """Convert wavenumber (cm⁻¹) to angular frequency (rad/s)."""
    return nu_cm1 * CM1_TO_RADS


def rads_to_cm1(omega_rads):
    """Convert angular frequency (rad/s) to wavenumber (cm⁻¹)."""
    return omega_rads / CM1_TO_RADS


def cm1_to_eV(nu_cm1):
    """Convert wavenumber (cm⁻¹) to energy (eV)."""
    return nu_cm1 * CM1_TO_EV


def eV_to_J(E_eV):
    """Convert energy (eV) to energy (J)."""
    return E_eV * EV_TO_J
