#!/bin/bash
#SBATCH --job-name=naph-eph3
#SBATCH --partition=CPU
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --array=1-176%20
#SBATCH --output=../logs/eph3_%A_%a.out
#SBATCH --error=../logs/eph3_%A_%a.err

set -e
ulimit -s unlimited
module load vasp/6.4.3-impi

BASE="/curie-home/linkh/remote-claude-2/naphthalene-transport-v2"
MAPFILE="$BASE/code/slurm/eph_dirs.txt"

DIRNAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MAPFILE")
DIR="$BASE/data/eph/$DIRNAME"

mkdir -p "$BASE/data/logs"
cd "$DIR"

echo "=== EPH $(date) Task $SLURM_ARRAY_TASK_ID/$SLURM_ARRAY_TASK_MAX ==="
echo "Dir: $DIRNAME"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun -np $SLURM_NTASKS vasp_std

echo "=== EPH Done $(date) ==="
