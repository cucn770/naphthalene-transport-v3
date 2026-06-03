# Naphthalene Transport v3

[![Python](https://img.shields.io/badge/Python-3.13-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Complete-brightgreen)]()

> **Holstein-Peierls charge transport calculation for naphthalene crystal**
>
> Transfer integrals + electron-phonon coupling + mobility from first-principles DFT
>
> Based on: Wang et al., *J. Chem. Phys.* **127**, 044506 (2007)

---

## Overview

This repository computes the charge carrier mobility in crystalline naphthalene (P2₁/a space group) using the **Holstein-Peierls model** with parameters derived from density functional theory (DFT).

### Pipeline (6 Stages)

```
Stage 1: TB Band Fitting → t_mn (7 transfer integrals, 2×2 model with P2₁/a symmetry)
Stage 2: DFPT Phonons    → ω_λ, e_λ (108 modes, 88 selected)
Stage 3: E-ph Coupling   → ∂ε₀/∂Q, ∂t/∂Q (176 displacement structures)
Stage 4: Parameter Assembly → G_λ (dimensionless coupling)
Stage 5: Mobility Prefactors → A₀, B (with E2 fix)
Stage 6: Time Integration → μ_α(T) (with E1 fix)
```

### v3 Corrections (E1-E4)

This version fixes 4 critical errors identified in the v2 pipeline:

| # | Error | Fix | Impact |
|---|-------|-----|--------|
| **E1** | Time integral missing factor 2 (even integrand) | `I₀, I_λ × 2.0` | μ × 2 |
| **E2** | Neighbor summation missing ±R paths | `A₀, B × 2.0` | μ × 2 |
| **E3** | H_AB(k) missing P2₁/a screw-axis symmetry | Use `build_tb_H_AB()` with `2·cos(k·τ)·exp(ik·R)` | t_ab: 60→30 meV |
| **E4** | Cross-sublattice sign degeneracy in finite difference | Continuity check vs equilibrium | 0/88 sign warnings |

**E1+E2 verified**: Pure paper-input test gives μ_c'(300K) = 0.964 vs v2's 0.241 = exactly 4.00×.

---

## Key Results

### Transfer Integrals

| Param | v3 (meV) | Paper (meV) |
|-------|----------|-------------|
| ε₀ | -388.3 | -835.0 |
| t_a | +37.03 | -23.00 |
| t_b | -28.86 | -42.00 |
| t_c | -5.30 | -3.00 |
| t_ac | -0.17 | -1.00 |
| t_ab | +30.49 | +22.00 |
| t_abc | +7.58 | -5.00 |
| **RMSE** | **26.15** | — |

### Mobility μ(300K)

| Direction | v3 (cm²/Vs) | Paper Calc. | Experiment (TOF) |
|-----------|-------------|-------------|------------------|
| a | **28.99** | — | ~1–5 |
| b | **10.62** | — | ~0.8–2 |
| c' | **1.62** | ~0.1 | ~0.4–1.0 |

Convergence: Δ = 0.0% (dt-halving test, all directions).

---

## Quick Start

### Prerequisites

- Python 3.10+ with NumPy, SciPy
- VASP 6.4.3 (for raw DFT data generation)
- Phonopy 4.1.0 (for DFPT post-processing)

### Installation

```bash
git clone https://github.com/cucn770/naphthalene-transport-v3.git
cd naphthalene-transport-v3
pip install numpy scipy phonopy pyyaml
```

### Data Setup

This repository includes computed results (`.npz` files) but NOT the raw VASP data (EIGENVAL, FORCE_CONSTANTS, etc.) due to size constraints.
See [`data/README.md`](data/README.md) for detailed setup instructions.

Quick setup with existing v2 project:

```bash
ln -s /path/to/naphthalene-transport-v2/data/crystal data/crystal
ln -s /path/to/naphthalene-transport-v2/data/phonon  data/phonon
ln -s /path/to/naphthalene-transport-v2/data/eph     data/eph
```

### Run the Pipeline

```bash
# Stage 1: TB Fitting
python3 -c "
import sys; sys.path.insert(0,'code')
from naph_transport.stage1_transfer import run_stage1
run_stage1(verbose=True)
"

# Stage 2: Phonon Post-processing (requires DFPT data)
# Stage 3: E-ph Coupling (requires 176 EIGENVAL files)
# Stage 4-6: Parameters + Mobility
python3 -c "
import sys; sys.path.insert(0,'code')
from naph_transport.stage4_assemble import run_stage4
from naph_transport.stage5_mobility_assemble import run_stage5
from naph_transport.stage6_mobility_integrate import run_stage6
run_stage4(verbose=True)
run_stage5(verbose=True)
run_stage6(T_min=10, T_max=300, T_step=10, verbose=True)
"

# Full Validation
python3 code/validate_pipeline.py
```

### Load Pre-computed Results

```python
import sys; sys.path.insert(0, 'code')
from naph_transport.io_utils import load_checkpoint
import numpy as np

# Transfer integrals
s1 = load_checkpoint('data/results/stage1_transfer_integrals.npz')
print(f"t_ab = {float(s1['t_ab'])*1000:.2f} meV, RMSE = {float(s1['rmse_eV'])*1000:.2f} meV")

# Mobility
s6 = load_checkpoint('data/results/stage6_mobility_results.npz')
T, mu_c = s6['T_range'], s6['mu_total_c_prime']
idx_300 = np.argmin(np.abs(T - 300))
print(f"mu_c'(300K) = {mu_c[idx_300]:.4f} cm²/Vs")
```

---

## Repository Structure

```
naphthalene-transport-v3/
├── code/                                  # Python package (v3.0.0)
│   ├── naph_transport/
│   │   ├── __init__.py                    # v3.0.0 with E1-E4 changelog
│   │   ├── constants.py                   # Physical constants + unit conversion
│   │   ├── crystal_utils.py               # P2₁/a crystal + TB basis functions
│   │   ├── io_utils.py                    # VASP/Phonopy I/O parsers
│   │   ├── validator.py                   # CP1-CP6 checkpoints
│   │   ├── stage1_transfer.py             # [E3] 2×2 TB fitting → t_mn
│   │   ├── stage2_phonon.py               # DFPT post-processing
│   │   ├── stage3_eph.py                  # [E4] E-ph coupling ∂t/∂Q, ∂ε₀/∂Q
│   │   ├── stage4_assemble.py             # G_λ parameter assembly
│   │   ├── stage5_mobility_assemble.py    # [E2] A₀, B prefactors
│   │   └── stage6_mobility_integrate.py   # [E1] μ(T) time integration
│   ├── validate_pipeline.py               # End-to-end validation
│   ├── slurm/                             # HPC job submission scripts
│   └── vasp_inputs/                       # VASP input templates
├── data/results/                          # Computed results (.npz files)
│   ├── stage1_transfer_integrals.npz
│   ├── stage2_phonon.npz
│   ├── stage3_eph_coupling.npz
│   ├── stage4_assembled.npz
│   ├── stage5_mobility_params.npz
│   └── stage6_mobility_results.npz
├── docs/                                  # Documentation
│   ├── v3_error_fixes.md
│   └── superpowers/
│       ├── specs/2026-06-03-naphthalene-v3-corrections-design.md
│       └── plans/2026-06-03-naphthalene-v3-corrections-plan.md
├── PIPELINE_COMPLETION_REPORT_V3.md       # Complete results report
├── TASK_RECORD.md                         # Full task history & file inventory
└── README.md                              # This file
```

---

## Validation (6/6 Passed)

```
CP1: TB Band Fitting ........ PASS (κ=2.4, RMSE=26.15 meV, E(k)=E(-k) ✓)
CP2: Phonon Spectrum ........ PASS (108 modes, no imaginary frequencies)
CP3: E-ph Coupling .......... PASS (88 modes, 0 sign-continuity warnings)
CP4: Parameter Assembly ..... PASS (G_total max=0.654, A_α OK)
CP5: Mobility Parameters .... PASS (A₀, B with E2 correction)
CP6: Mobility Results ....... PASS (Δ=0.0% convergence, all directions)
```

---

## References

1. Wang, L. J. et al. "Electron-vibration coupling in the charge transport of naphthalene." *J. Chem. Phys.* **127**, 044506 (2007). [DOI: 10.1063/1.2751191](https://doi.org/10.1063/1.2751191)

2. Senthilkumar, K. et al. "Charge transport in organic semiconductors." *Phys. Rev. Lett.* **96**, 086601 (2006). [DOI: 10.1103/PhysRevLett.96.086601](https://doi.org/10.1103/PhysRevLett.96.086601)

3. Valeev, E. F. et al. "Effect of electronic polarization on charge-transport parameters in molecular organic semiconductors." *J. Am. Chem. Soc.* **128**, 9882 (2006). [DOI: 10.1021/ja061827h](https://doi.org/10.1021/ja061827h)

4. Karl, N. "Organic Semiconductors." In *Festkörperprobleme XIV*, 261 (1974).

5. Warta, W. & Karl, N. "Hot holes in naphthalene." *Phys. Rev. B* **32**, 1172 (1985).

---

## License

MIT License — see [LICENSE](LICENSE) file for details.

## Citation

If you use this code in your research, please cite both this repository and the original Wang et al. (2007) paper.
