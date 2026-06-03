"""Crystal structure utilities for naphthalene P2_1/a.

Provides lattice vectors, k-point conversion (fractional→Cartesian),
and tight-binding basis function construction for the 2×2 model.
"""
import numpy as np
from dataclasses import dataclass
from typing import Tuple, Dict


@dataclass
class CrystalParams:
    """Naphthalene crystal parameters (CSD experimental, Ponomarev 1976)."""
    a: float = 8.098      # Å
    b: float = 5.953      # Å
    c: float = 8.652      # Å
    beta_deg: float = 124.4
    space_group: str = "P2_1/a"
    n_molecules_per_cell: int = 2
    n_atoms_per_molecule: int = 18
    n_atoms_total: int = 36

    @property
    def beta_rad(self) -> float:
        return np.radians(self.beta_deg)

    def lattice_vectors(self) -> np.ndarray:
        """Return (3,3) real-space lattice vectors in Å (Cartesian)."""
        return np.array([
            [self.a, 0.0, 0.0],
            [0.0, self.b, 0.0],
            [self.c * np.cos(self.beta_rad), 0.0, self.c * np.sin(self.beta_rad)],
        ])

    def reciprocal_vectors(self) -> np.ndarray:
        """Return (3,3) reciprocal lattice vectors B = [b1; b2; b3] in Å⁻¹ (Cartesian)."""
        lat = self.lattice_vectors()
        volume = np.dot(lat[0], np.cross(lat[1], lat[2]))
        recip = np.zeros((3, 3))
        recip[0] = 2.0 * np.pi * np.cross(lat[1], lat[2]) / volume
        recip[1] = 2.0 * np.pi * np.cross(lat[2], lat[0]) / volume
        recip[2] = 2.0 * np.pi * np.cross(lat[0], lat[1]) / volume
        return recip

    def k_frac_to_cart(self, k_frac: np.ndarray) -> np.ndarray:
        """Convert fractional k-point coordinates to Cartesian (Å⁻¹).

        Args:
            k_frac: (N, 3) or (3,) fractional coordinates

        Returns:
            k_cart: (N, 3) or (3,) Cartesian coordinates in Å⁻¹
        """
        B = self.reciprocal_vectors()  # (3,3)
        return k_frac @ B

    def R_frac_to_cart(self, R_frac: np.ndarray) -> np.ndarray:
        """Convert fractional translation vector to Cartesian (Å).

        Args:
            R_frac: (3,) fractional coordinates

        Returns:
            R_cart: (3,) Cartesian coordinates in Å
        """
        lat = self.lattice_vectors()
        return R_frac @ lat


def build_tb_H_AA(k_cart: np.ndarray, R_cart_same: np.ndarray, t_same: np.ndarray) -> np.ndarray:
    """Build H_AA(k) = 2 * Σ t_n * cos(k·R_n) for same-sublattice hopping.

    Args:
        k_cart: (N_k, 3) k-point coordinates in Cartesian (Å⁻¹)
        R_cart_same: (N_same, 3) translation vectors in Cartesian (Å)
        t_same: (N_same,) transfer integrals in eV

    Returns:
        H_AA: (N_k,) same-sublattice contribution in eV
    """
    H_AA = np.zeros(k_cart.shape[0])
    for i in range(len(t_same)):
        k_dot_R = np.dot(k_cart, R_cart_same[i])
        H_AA += 2.0 * t_same[i] * np.cos(k_dot_R)
    return H_AA


def build_tb_H_AB(k_cart: np.ndarray, R_cart_cross: np.ndarray,
                   t_cross: np.ndarray, tau_AB: np.ndarray) -> np.ndarray:
    """Build H_AB(k) for cross-sublattice hopping (complex phase factors).

    For P2_1/a with molecules A and B at relative displacement tau_AB,
    the cross-sublattice hopping for each translation R has contributions
    from the 2 screw-axis-equivalent positions of B relative to A:

        H_AB(k) = Σ_n t_n · [exp(i·k·(R_n+τ_AB)) + exp(i·k·(R_n-τ_AB))]

    which simplifies to 2·t_n·cos(k·τ_AB)·exp(i·k·R_n).

    Args:
        k_cart: (N_k, 3) k-point coordinates in Cartesian (Å⁻¹)
        R_cart_cross: (N_cross, 3) representative translation vectors in Cartesian (Å)
        t_cross: (N_cross,) transfer integrals in eV
        tau_AB: (3,) relative displacement A→B within unit cell in Cartesian (Å)

    Returns:
        H_AB_complex: (N_k,) complex cross-sublattice contribution in eV
    """
    H_AB = np.zeros(k_cart.shape[0], dtype=complex)

    for i in range(len(t_cross)):
        R = R_cart_cross[i]
        k_dot_R = np.dot(k_cart, R)
        k_dot_tau = np.dot(k_cart, tau_AB)

        # Sum over ±tau_AB for P2_1/a screw axis
        # t_n * [exp(i·k·(R+τ)) + exp(i·k·(R-τ))] = 2·t_n·cos(k·τ)·exp(i·k·R)
        H_AB += t_cross[i] * 2.0 * np.cos(k_dot_tau) * np.exp(1j * k_dot_R)

    return H_AB


def tb_eigenvalues(k_cart: np.ndarray, epsilon_0: float,
                   t_same: Dict[str, float], R_cart_same: np.ndarray,
                   t_cross: Dict[str, float], R_cart_cross: np.ndarray,
                   tau_AB: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Compute the two TB eigenvalues E±(k) for the 2×2 model.

    H(k) = [[ε₀ + H_AA,   H_AB      ],
            [H_AB*,        ε₀ + H_AA ]]

    E±(k) = ε₀ + H_AA(k) ± |H_AB(k)|

    Args:
        k_cart: (N_k, 3) k-points in Cartesian (Å⁻¹)
        epsilon_0: on-site energy (eV)
        t_same: dict of same-sublattice t values {'a': t_a, 'b': t_b, 'c': t_c}
        R_cart_same: (3, 3) same-sublattice R vectors
        t_cross: dict of cross-sublattice t values {'ac': t_ac, 'ab': t_ab, 'abc': t_abc}
        R_cart_cross: (3, 3) cross-sublattice R vectors
        tau_AB: (3,) sublattice displacement in Cartesian (Å)

    Returns:
        E_plus: (N_k,) upper band (HOMO) energies in eV
        E_minus: (N_k,) lower band (HOMO-1) energies in eV
    """
    t_same_arr = np.array([t_same['a'], t_same['b'], t_same['c']])
    t_cross_arr = np.array([t_cross['ac'], t_cross['ab'], t_cross['abc']])

    H_AA = build_tb_H_AA(k_cart, R_cart_same, t_same_arr)
    H_AB = build_tb_H_AB(k_cart, R_cart_cross, t_cross_arr, tau_AB)

    H_AB_abs = np.abs(H_AB)

    E_plus = epsilon_0 + H_AA + H_AB_abs
    E_minus = epsilon_0 + H_AA - H_AB_abs

    return E_plus, E_minus


def compute_sublattice_displacement(positions_frac: np.ndarray,
                                     lattice: np.ndarray) -> np.ndarray:
    """Compute the relative displacement vector tau_AB between sublattices.

    Takes the center-of-mass of molecule A and molecule B to find tau_AB.

    Args:
        positions_frac: (36, 3) fractional coordinates
        lattice: (3, 3) lattice vectors in Å

    Returns:
        tau_AB: (3,) Cartesian displacement A→B in Å
    """
    from .constants import MOL_A_ATOMS, MOL_B_ATOMS, MASS_C, MASS_H

    # Build mass array: first 10 C of mol A, then 8 H of mol A, then 10 C of mol B, then 8 H of mol B
    masses = np.array([MASS_C]*10 + [MASS_H]*8 + [MASS_C]*10 + [MASS_H]*8)

    # Convert fractional to Cartesian
    pos_cart = positions_frac @ lattice

    # Compute center of mass for each sublattice
    mol_a_mass = np.sum(masses[MOL_A_ATOMS])
    mol_b_mass = np.sum(masses[MOL_B_ATOMS])

    com_a = np.sum(pos_cart[MOL_A_ATOMS] * np.array(masses)[MOL_A_ATOMS][:, np.newaxis], axis=0) / mol_a_mass
    com_b = np.sum(pos_cart[MOL_B_ATOMS] * np.array(masses)[MOL_B_ATOMS][:, np.newaxis], axis=0) / mol_b_mass

    return com_b - com_a
