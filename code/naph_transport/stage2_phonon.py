"""Stage 2: Post-process Phonopy DFPT results.

Reads FORCE_CONSTANTS from Phonopy, extracts frequencies at Gamma,
selects modes for Stage 3 e-ph coupling (exclude acoustic + >2800 cm⁻¹),
generates displacement structures using Phonopy MODULATION.
"""
import numpy as np
import subprocess
import os
from pathlib import Path
from typing import Dict, Tuple
import yaml

from .io_utils import save_checkpoint
from .validator import cp2_check_phonon


def read_phonon_frequencies(band_yaml_path: str = 'data/phonon/band.yaml'
                            ) -> np.ndarray:
    """Read phonon frequencies at Gamma from Phonopy band.yaml.

    Args:
        band_yaml_path: path to band.yaml

    Returns:
        frequencies_cm1: (N_modes,) frequencies at Gamma (cm⁻¹)
    """
    with open(band_yaml_path, 'r') as f:
        data = yaml.safe_load(f)

    for qpt in data['phonon']:
        q_pos = np.array(qpt['q-position'])
        if np.allclose(q_pos, 0.0, atol=1e-6):
            return np.array([band['frequency'] for band in qpt['band']])

    raise ValueError("Gamma point not found in band.yaml")


def select_modes_for_eph(frequencies_cm1: np.ndarray,
                         freq_min: float = 30.0,
                         freq_max: float = 2800.0
                         ) -> Tuple[np.ndarray, np.ndarray]:
    """Select optical modes for e-ph coupling calculation.

    Excludes: 3 acoustic modes (ω ≈ 0) and pure C-H stretches (>2800 cm⁻¹).

    Args:
        frequencies_cm1: (108,) all mode frequencies
        freq_min: minimum frequency for selection (cm⁻¹)
        freq_max: maximum frequency for selection (cm⁻¹)

    Returns:
        selected_indices: 0-indexed mode indices
        selected_frequencies: corresponding frequencies
    """
    # Sort |ω| to identify acoustic modes
    sorted_idx = np.argsort(np.abs(frequencies_cm1))
    acoustic_idx = sorted_idx[:3]

    # Select all other modes within frequency range
    mask = np.ones(len(frequencies_cm1), dtype=bool)
    mask[acoustic_idx] = False
    mask = mask & (np.abs(frequencies_cm1) >= freq_min)
    mask = mask & (np.abs(frequencies_cm1) <= freq_max)

    selected_indices = np.where(mask)[0]
    selected_frequencies = frequencies_cm1[selected_indices]

    return selected_indices, selected_frequencies


def generate_displacements(mode_index: int, amplitude: float = 0.01,
                           phonon_dir: str = 'data/phonon',
                           eph_dir: str = 'data/eph') -> Tuple[str, str]:
    """Generate ±ΔQ displacement structures for a given mode using Phonopy.

    Uses Phonopy's MODULATION tag to ensure correct mass-weighted displacements.

    Args:
        mode_index: 0-indexed phonon mode index
        amplitude: displacement amplitude in √amu·Å
        phonon_dir: directory with FORCE_CONSTANTS and POSCAR
        eph_dir: output directory

    Returns:
        dir_plus, dir_minus: paths to displacement directories
    """
    mode_tag = f"mode_{mode_index+1:04d}"
    dir_plus = os.path.join(eph_dir, f"eph_{mode_tag}_plus")
    dir_minus = os.path.join(eph_dir, f"eph_{mode_tag}_minus")

    for sign, outdir in [('+', dir_plus), ('-', dir_minus)]:
        os.makedirs(outdir, exist_ok=True)

        mod_amplitude = amplitude if sign == '+' else -amplitude

        # Use phonopy command to generate modulation
        cmd = [
            'phonopy', '--dim', '1', '1', '1',
            '-c', os.path.join(phonon_dir, 'POSCAR'),
            '--readfc',
            '--modulation', f'{mode_index+1} {mode_index+1} {mod_amplitude}',
        ]
        result = subprocess.run(cmd, cwd=phonon_dir, capture_output=True, text=True)

        # Phonopy writes MPOSCAR in the current directory
        mposcar = os.path.join(phonon_dir, 'MPOSCAR')
        if os.path.exists(mposcar):
            import shutil
            shutil.move(mposcar, os.path.join(outdir, 'POSCAR'))
        else:
            print(f"WARNING: MPOSCAR not generated for mode {mode_index+1} {sign}")
            print(f"stdout: {result.stdout}")
            print(f"stderr: {result.stderr}")

    return dir_plus, dir_minus


def run_stage2(phonon_dir: str = 'data/phonon',
               verbose: bool = True) -> Dict:
    """Run Stage 2: post-process DFPT phonons and select modes.

    Args:
        phonon_dir: directory with DFPT output (FORCE_CONSTANTS, band.yaml)
        verbose: print progress

    Returns:
        dict with frequencies, selected mode indices, and counts
    """
    # Generate band.yaml if not present
    band_yaml = os.path.join(phonon_dir, 'band.yaml')
    if not os.path.exists(band_yaml):
        subprocess.run([
            'phonopy', '--dim', '1', '1', '1',
            '-c', os.path.join(phonon_dir, 'POSCAR'),
            '--readfc', '--band', '0 0 0  0.5 0 0',
            '--band-points', '51'
        ], cwd=phonon_dir, check=True)

    frequencies = read_phonon_frequencies(band_yaml)

    if verbose:
        print(f"Stage 2: Phonon post-processing")
        print(f"  Total modes: {len(frequencies)}")
        ac = np.sort(np.abs(frequencies[:3]))
        print(f"  Acoustic (|ω|): {ac} cm⁻¹")
        print(f"  Optical range: {np.min(np.abs(frequencies[3:])):.1f} to "
              f"{np.max(frequencies):.1f} cm⁻¹")

    # Validate
    results = cp2_check_phonon(frequencies)

    # Select modes for Stage 3
    selected_idx, selected_freq = select_modes_for_eph(frequencies)

    if verbose:
        print(f"  Selected for e-ph: {len(selected_idx)} modes "
              f"({np.min(selected_freq):.1f}–{np.max(selected_freq):.1f} cm⁻¹)")
        print(f"  CP2: {'ALL PASS' if all(results.values()) else 'see details above'}")

    # Save
    save_checkpoint({
        'frequencies_cm1': frequencies,
        'selected_indices': selected_idx,
        'selected_frequencies': selected_freq,
    }, 'data/results/stage2_phonon.npz')

    return {
        'frequencies_cm1': frequencies,
        'selected_indices': selected_idx,
        'selected_frequencies': selected_freq,
        'n_selected': len(selected_idx),
    }
