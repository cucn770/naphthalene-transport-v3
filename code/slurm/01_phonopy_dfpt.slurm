#!/bin/bash
#SBATCH --job-name=naph-dfpt
#SBATCH --partition=CPU
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --output=../logs/dfpt_%j.out
#SBATCH --error=../logs/dfpt_%j.err

# DFPT phonon for naphthalene (36 atoms, P2_1/a)
# IBRION=8, LEPSILON=.TRUE.
# 36 atoms × 3 DOF = 108 perturbations
# 32 cores on one node, walltime 24h

set -e
ulimit -s unlimited
module load vasp/6.4.3-impi

WORKDIR="/curie-home/linkh/remote-claude-2/naphthalene-transport-v2/data/phonon"
mkdir -p "$WORKDIR/../logs"
cd "$WORKDIR"

echo "=== DFPT Start: $(date) ==="
echo "Host: $(hostname)"
echo "Working: $(pwd)"
echo "Cores: $SLURM_NTASKS"

# Run VASP DFPT
mpirun -np $SLURM_NTASKS vasp_std

echo "=== DFPT End: $(date) ==="

# Post-process with Phonopy
cat > band.conf << 'EOF'
DIM = 1 1 1
BAND = 0.0 0.0 0.0  0.5 0.0 0.0  0.0 0.5 0.0  0.0 0.0 0.5  0.0 0.0 0.0
BAND_LABELS = Γ X Γ Y Γ Z Γ
BAND_POINTS = 101
FULL_FORCE_CONSTANTS = .TRUE.
FC_SYMMETRY = .TRUE.
EOF

cat > mesh.conf << 'EOF'
DIM = 1 1 1
MP = 4 4 4
FULL_FORCE_CONSTANTS = .TRUE.
FC_SYMMETRY = .TRUE.
EOF

# Generate FORCE_CONSTANTS from vasprun.xml
phonopy --fc vasprun.xml

echo "=== Phonopy band structure ==="
phonopy -p band.conf 2>&1 | head -30

echo "=== DFPT workflow complete: $(date) ==="
