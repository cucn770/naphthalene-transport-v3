#!/bin/bash
# Master submit script for naphthalene transport v2 pipeline
# Usage: bash master_submit.sh [stage]

set -e

BASE="/curie-home/linkh/remote-claude-2/naphthalene-transport-v2"
STAGE="${1:-all}"

echo "=== Naphthalene Transport v2 Pipeline ==="
echo "Working: $BASE"
echo "Stage: $STAGE"

submit_stage1() {
    echo "--- Stage 1: TB Band Fitting ---"
    cd "$BASE"
    python3 -c "
import sys; sys.path.insert(0, 'code')
from naph_transport.stage1_transfer import run_stage1
t = run_stage1()
" && echo "Stage 1 complete." || echo "Stage 1 FAILED."
}

submit_stage2() {
    echo "--- Stage 2: Phonopy DFPT ---"
    cd "$BASE/data/phonon"
    JOBID=$(sbatch --parsable ../../code/slurm/01_phonopy_dfpt.slurm)
    echo "DFPT submitted: Job $JOBID"
    echo "Monitor: squeue -j $JOBID"
}

submit_stage2_post() {
    echo "--- Stage 2 Post: Phonopy Analysis ---"
    cd "$BASE"
    python3 -c "
import sys; sys.path.insert(0, 'code')
from naph_transport.stage2_phonon import run_stage2
result = run_stage2()
print(f'Selected {result[\"n_selected\"]} modes for Stage 3')
" && echo "Stage 2 post complete." || echo "Stage 2 post FAILED."
}

submit_stage3_setup() {
    echo "--- Stage 3 Setup: Generate Displacement Directories ---"
    cd "$BASE"
    python3 << 'PYEOF'
import sys, os, shutil
sys.path.insert(0, 'code')
from naph_transport.io_utils import load_checkpoint
from naph_transport.stage2_phonon import generate_displacements

s2 = load_checkpoint('data/results/stage2_phonon.npz')
selected_idx = [int(x) for x in s2['selected_indices']]
n_modes = len(selected_idx)
print(f'Generating {n_modes} modes × 2 = {2*n_modes} displacement dirs...')

for i, mode_idx in enumerate(selected_idx):
    dir_plus, dir_minus = generate_displacements(
        int(mode_idx), amplitude=0.01,
        phonon_dir='data/phonon', eph_dir='data/eph')

    # Setup VASP inputs for each displacement dir
    for d in [dir_plus, dir_minus]:
        if os.path.exists(os.path.join(d, 'POSCAR')):
            shutil.copy('code/vasp_inputs/INCAR_eph', os.path.join(d, 'INCAR'))
            shutil.copy('code/vasp_inputs/KPOINTS_444', os.path.join(d, 'KPOINTS'))
            shutil.copy('data/crystal/POTCAR', os.path.join(d, 'POTCAR'))

print(f'All {2*n_modes} displacement directories prepared.')
print(f'Update 02_eph_array.slurm: --array=1-{2*n_modes}%20')
PYEOF
    echo "Stage 3 setup complete."
}

submit_stage3() {
    echo "--- Stage 3: E-ph Coupling Array ---"
    cd "$BASE/data/eph"
    # Count actual directories to set correct array size
    N_DIRS=$(ls -d eph_mode_* 2>/dev/null | wc -l)
    if [ "$N_DIRS" -eq 0 ]; then
        echo "ERROR: No displacement directories found. Run stage3_setup first."
        exit 1
    fi
    echo "Found $N_DIRS displacement directories."
    # Submit with correct array range
    JOBID=$(sbatch --parsable --array="1-${N_DIRS}%20" ../../code/slurm/02_eph_array.slurm)
    echo "E-ph array submitted: Job $JOBID"
    echo "Monitor: squeue -j $JOBID"
}

submit_stage3_post() {
    echo "--- Stage 3 Post: E-ph Analysis ---"
    cd "$BASE"
    python3 << 'PYEOF'
import sys, os, glob
sys.path.insert(0, 'code')
from naph_transport.io_utils import load_checkpoint, save_checkpoint
from naph_transport.stage3_eph import run_stage3
import numpy as np

# Load equilibrium TB parameters
s1 = load_checkpoint('data/results/stage1_transfer_integrals.npz')
t_eq = {
    'epsilon_0': float(s1['epsilon_0']),
    'a': float(s1['t_a']), 'b': float(s1['t_b']), 'c': float(s1['t_c']),
    'ac': float(s1['t_ac']), 'ab': float(s1['t_ab']), 'abc': float(s1['t_abc']),
}
tau_AB = s1['tau_AB']

# Find all displacement directories
eph_dirs_plus = sorted(glob.glob('data/eph/eph_mode_*_plus'))
eph_dirs_minus = sorted(glob.glob('data/eph/eph_mode_*_minus'))

# Pair them (plus, minus alternating)
eph_dirs = []
for p, m in zip(eph_dirs_plus, eph_dirs_minus):
    eph_dirs.extend([p, m])

print(f'Processing {len(eph_dirs)//2} modes from {len(eph_dirs)} directories...')

result = run_stage3(eph_dirs, t_eq, tau_AB, delta_Q=0.01, verbose=True)
save_checkpoint({
    'd_epsilon0_dQ': result['d_epsilon0_dQ'],
    'dt_dQ': result['dt_dQ'],
}, 'data/results/stage3_eph_coupling.npz')
print('Stage 3 complete.')
PYEOF
    echo "Stage 3 post complete."
}

submit_stage4() {
    echo "--- Stage 4: Parameter Assembly ---"
    cd "$BASE"
    python3 -c "
import sys; sys.path.insert(0, 'code')
from naph_transport.stage4_assemble import run_stage4
result = run_stage4()
" && echo "Stage 4 complete." || echo "Stage 4 FAILED."
}

case "$STAGE" in
    all)
        echo "=== Full Pipeline ==="
        submit_stage1
        echo ""
        echo "=== Next: submit DFPT (Stage 2) manually ==="
        echo "  bash $0 stage2"
        ;;
    1) submit_stage1 ;;
    2) submit_stage2 ;;
    2post) submit_stage2_post ;;
    3setup) submit_stage3_setup ;;
    3) submit_stage3 ;;
    3post) submit_stage3_post ;;
    4) submit_stage4 ;;
    pipeline)
        submit_stage1
        submit_stage2
        echo ""
        echo "=== After DFPT completes: ==="
        echo "  bash $0 2post"
        echo "  bash $0 3setup"
        echo "  bash $0 3"
        echo "  bash $0 3post"
        echo "  bash $0 4"
        ;;
    *)
        echo "Usage: bash master_submit.sh [1|2|2post|3setup|3|3post|4|all|pipeline]"
        exit 1
        ;;
esac
