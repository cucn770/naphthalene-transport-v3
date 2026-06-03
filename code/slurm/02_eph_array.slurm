#!/bin/bash
#SBATCH --job-name=naph-eph
#SBATCH --partition=CPU
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --array=1-160%20
#SBATCH --output=../logs/eph_%A_%a.out
#SBATCH --error=../logs/eph_%A_%a.err

# E-ph coupling: single-point band calculation for displaced structure
# Array job: each task runs one displacement direction
# Uses WAVECAR from equilibrium structure for 3-5 step SCF

set -e
ulimit -s unlimited
module load vasp/6.4.3-impi

BASE="/curie-home/linkh/remote-claude-2/naphthalene-transport-v2"
CRYSTAL="$BASE/data/crystal"
EPH_BASE="$BASE/data/eph"

# Determine mode number and sign from array index
# Array indices 1..N: even=plus, odd=minus
ARRAY_IDX=$((SLURM_ARRAY_TASK_ID - 1))
MODE_IDX=$((ARRAY_IDX / 2 + 1))
SIGN_IDX=$((ARRAY_IDX % 2))

MODE_TAG=$(printf "mode_%04d" $MODE_IDX)
if [ $SIGN_IDX -eq 0 ]; then
    DIR="$EPH_BASE/eph_${MODE_TAG}_plus"
else
    DIR="$EPH_BASE/eph_${MODE_TAG}_minus"
fi

mkdir -p "$BASE/data/logs"
cd "$DIR"

echo "=== EPH $(date) ==="
echo "Array ID: $SLURM_ARRAY_TASK_ID, Mode: $MODE_TAG, Sign: $([ $SIGN_IDX -eq 0 ] && echo '+' || echo '-')"
echo "Working: $DIR"

# Link input files
ln -sf "$CRYSTAL/POTCAR" POTCAR 2>/dev/null || cp "$CRYSTAL/POTCAR" .
cp "$BASE/code/vasp_inputs/KPOINTS_444" KPOINTS

# Copy WAVECAR and CHGCAR from equilibrium crystal for SCF acceleration
cp "$CRYSTAL/WAVECAR" .
cp "$CRYSTAL/CHGCAR" .

# Run VASP (4 MPI × 2 OMP = 8 cores)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun -np $SLURM_NTASKS vasp_std

# Quick validation: check SCF convergence
if grep -q "reached required accuracy" OUTCAR; then
    echo "SCF converged successfully"
else
    echo "WARNING: SCF may not have converged — check OUTCAR"
fi

echo "=== EPH Done $(date) ==="
