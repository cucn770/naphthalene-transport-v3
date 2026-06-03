# Data Directory

This directory contains symbolic links to raw VASP calculation data from the v2 project.
These files are NOT included in the GitHub repository.

## Required Data for Full Pipeline Execution

To reproduce the full calculation, the following data must be available locally:

### Crystal Structure and DFT Bands
- `data/crystal/CONTCAR` — Optimized naphthalene crystal structure (P2₁/a)
- `data/crystal/band/EIGENVAL` — VASP band structure (36 k-points, 4×4×4 MP grid)
- `data/crystal/POTCAR` — PAW pseudopotentials (C, H)
- `data/crystal/WAVECAR` — Wavefunction (optional, for SCF acceleration)
- `data/crystal/CHGCAR` — Charge density (optional, for SCF acceleration)

### DFPT Phonon Data
- `data/phonon/FORCE_CONSTANTS` — DFPT force constants (262 KB, Phonopy format)
- `data/phonon/POSCAR` — Structure file used for DFPT
- `data/phonon/vasprun.xml` — VASP DFPT output

### Electron-Phonon Displacement Data
- `data/eph/eph_mode_XXXX_plus/EIGENVAL` — 88 files (36 k-point bands at +ΔQ)
- `data/eph/eph_mode_XXXX_minus/EIGENVAL` — 88 files (36 k-point bands at -ΔQ)
- `data/eph/eph_mode_XXXX_plus/OUTCAR` — 88 files (forces at +ΔQ)
- `data/eph/eph_mode_XXXX_minus/OUTCAR` — 88 files (forces at -ΔQ)

## Quick Setup (if you have v2 project locally)

```bash
ln -s /path/to/naphthalene-transport-v2/data/crystal data/crystal
ln -s /path/to/naphthalene-transport-v2/data/phonon  data/phonon
ln -s /path/to/naphthalene-transport-v2/data/eph     data/eph
ln -s /path/to/naphthalene-transport-v2/original_docs original_docs
```

## References

- Original VASP calculations: naphthalene-transport-v2 project
- DFPT: VASP 6.4.3, IBRION=8, LEPSILON=.TRUE., 32 CPU cores
- e-ph: VASP 6.4.3, IBRION=-1, NSW=0, 176 jobs (88 modes × ±ΔQ)
- Computational details: see `docs/CALCULATION_WORKFLOW.md` in v2 project
