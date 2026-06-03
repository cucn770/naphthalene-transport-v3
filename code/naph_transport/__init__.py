"""naph_transport: Naphthalene charge transport parameter calculation.

v3.0.0: Fixed 4 critical errors in the Holstein-Peierls mobility pipeline
  - E1: Time integral even-function factor of 2 (Stage 6)
  - E2: Neighbor summation +/-R factor of 2 (Stage 5)
  - E3: H_AB(k) P2_1/a screw-axis symmetry (Stage 1)
  - E4: Cross-sublattice sign continuity check (Stage 3)
"""
__version__ = "3.0.0"
