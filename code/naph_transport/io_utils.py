"""I/O utilities for parsing VASP, Phonopy, and CONTCAR files."""
import numpy as np
from pathlib import Path
from typing import Tuple, Dict


def parse_vasp_eigenval(eigenval_path: str) -> Tuple[np.ndarray, np.ndarray, int]:
    """Parse VASP EIGENVAL file to extract k-points and band energies.

    Args:
        eigenval_path: path to EIGENVAL file

    Returns:
        k_points: (N_k, 3) k-point coordinates (fractional reciprocal)
        bands: (N_k, N_bands) band energies (eV)
        n_electrons: number of electrons
    """
    with open(eigenval_path, 'r') as f:
        lines = f.readlines()

    # Line 6: #electrons, #kpoints, #bands
    header_parts = lines[5].strip().split()
    n_electrons = int(header_parts[0])
    n_kpoints = int(header_parts[1])
    n_bands = int(header_parts[2])

    k_points = np.zeros((n_kpoints, 3))
    bands = np.zeros((n_kpoints, n_bands))

    idx = 7
    for ik in range(n_kpoints):
        while idx < len(lines) and lines[idx].strip() == '':
            idx += 1
        if idx >= len(lines):
            break

        parts = lines[idx].strip().split()
        k_points[ik] = [float(x) for x in parts[:3]]
        idx += 1

        for ib in range(n_bands):
            if idx >= len(lines):
                break
            parts = lines[idx].strip().split()
            bands[ik, ib] = float(parts[1])
            idx += 1

    return k_points, bands, n_electrons


def parse_vasp_contcar(contcar_path: str) -> Tuple[np.ndarray, np.ndarray, list]:
    """Parse VASP CONTCAR to get lattice vectors and atomic positions.

    Args:
        contcar_path: path to CONTCAR file

    Returns:
        lattice: (3, 3) lattice vectors (Å)
        positions: (N_atoms, 3) atomic positions (direct/fractional coordinates)
        atom_symbols: list of atom symbols
    """
    with open(contcar_path, 'r') as f:
        lines = f.readlines()

    scale = float(lines[1].strip())

    lattice = np.array([
        [float(x) for x in lines[2].strip().split()],
        [float(x) for x in lines[3].strip().split()],
        [float(x) for x in lines[4].strip().split()],
    ]) * scale

    line6 = lines[5].strip()
    line7 = lines[6].strip()
    try:
        n_atoms_per_type = [int(x) for x in line7.split()]
        atom_symbols_list = line6.split()
    except ValueError:
        n_atoms_per_type = [int(x) for x in line6.split()]
        atom_symbols_list = line7.split()

    n_atoms = sum(n_atoms_per_type)

    positions = np.zeros((n_atoms, 3))
    start_line = 8
    for i in range(n_atoms):
        parts = lines[start_line + i].strip().split()
        positions[i] = [float(x) for x in parts[:3]]

    return lattice, positions, atom_symbols_list


def parse_vasp_outcar_forces(outcar_path: str) -> np.ndarray:
    """Parse final forces from VASP OUTCAR.

    Args:
        outcar_path: path to OUTCAR file

    Returns:
        forces: (N_atoms, 3) forces (eV/Å)
    """
    forces = []
    in_forces = False
    skip_separator = False
    force_block = []

    with open(outcar_path, 'r') as f:
        for line in f:
            if 'TOTAL-FORCE' in line:
                in_forces = True
                skip_separator = True
                force_block = []
                continue
            if in_forces:
                if skip_separator and '--------' in line:
                    skip_separator = False
                    continue
                if '--------' in line or line.strip() == '':
                    if len(force_block) > 0:
                        forces = force_block
                    in_forces = False
                    continue
                parts = line.strip().split()
                if len(parts) >= 6:
                    force_block.append([float(x) for x in parts[3:6]])

    return np.array(forces)


def save_checkpoint(data: dict, filepath: str) -> None:
    """Save checkpoint data as .npz file.

    Args:
        data: dictionary of numpy arrays or scalars
        filepath: path to save file
    """
    path = Path(filepath)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, **data)
    print(f"Checkpoint saved: {path}")


def load_checkpoint(filepath: str) -> dict:
    """Load checkpoint data from .npz file.

    Args:
        filepath: path to checkpoint file

    Returns:
        data: dictionary of numpy arrays
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Checkpoint not found: {filepath}")
    return dict(np.load(path))
