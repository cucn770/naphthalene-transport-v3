#!/usr/bin/env python3
"""Fix POSCAR: remove complex j from coordinates (eigenvectors must be .real at Gamma)."""
import sys; sys.path.insert(0, '.')
import numpy as np, os, shutil, glob
from collections import Counter
import phonopy
from phonopy.interface.vasp import read_vasp

unitcell = read_vasp("data/phonon/POSCAR")
phonon = phonopy.load(supercell_matrix=[[1,0,0],[0,1,0],[0,0,1]],
    primitive_matrix='P', unitcell=unitcell,
    force_constants_filename="data/phonon/FORCE_CONSTANTS", symprec=1e-5)

data = np.load('data/phonon/phonon_dfpt_results.npz')
fa = np.abs(data['frequencies_cm1'])
si = np.argsort(fa)
mask = np.ones(108, bool)
mask[si[:3]] = False  # exclude acoustic
mask &= (fa >= 30) & (fa <= 2800)
sel = np.where(mask)[0]

dm = phonon.dynamical_matrix; dm.run([0.,0.,0.])
_, ev = np.linalg.eigh(dm.dynamical_matrix)
ev = np.real(ev)  # CRITICAL FIX: Gamma eigenvectors are real

la = unitcell.cell; ila = np.linalg.inv(la)
pf = unitcell.scaled_positions; ms = unitcell.masses.reshape(-1,1)
sy = unitcell.symbols; sc = Counter(sy); us = sorted(sc, key=lambda x: sy.index(x))

cnt = 0
for mi in sel:
    mt = f"mode_{int(mi)+1:04d}"
    evm = ev[:, mi].reshape(36, 3)
    for ss, sd in [('plus', 0.01), ('minus', -0.01)]:
        dn = f"data/eph/eph_{mt}_{ss}"; os.makedirs(dn, exist_ok=True)
        for old in glob.glob(f"{dn}/*"): os.remove(old)
        pc = pf @ la
        pn = (pc + sd * evm / np.sqrt(ms)) @ ila
        with open(f"{dn}/POSCAR", 'w') as f:
            f.write(f"Naphthalene eph {mt} {ss}\n1.0\n")
            for v in la: f.write(f"  {v[0]:16.12f}  {v[1]:16.12f}  {v[2]:16.12f}\n")
            f.write(" ".join(us)+"\n"+ " ".join(str(sc[s]) for s in us)+"\n")
            f.write("Direct\n")
            for p in pn: f.write(f"  {float(p[0]):16.12f}  {float(p[1]):16.12f}  {float(p[2]):16.12f}\n")
        shutil.copy('code/vasp_inputs/INCAR_eph', f"{dn}/INCAR")
        shutil.copy('code/vasp_inputs/KPOINTS_444', f"{dn}/KPOINTS")
        shutil.copy('data/crystal/POTCAR', f"{dn}/POTCAR")
        cnt += 1
print(f"Done: {cnt} POSCARs regenerated")
