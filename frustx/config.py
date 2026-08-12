"""Tunable defaults, kept free of heavy imports.

This module exists so `frustx --help` does not have to load PyRosetta. The CLI needs
these values to build its argument parser; importing them from decoys.py or
frustration.py would drag in PyRosetta (~2.5 s) just to print usage text.

Nothing here may import pyrosetta, numpy or Bio.
"""

# Ca-Ca contact cutoff in Angstrom.
#   "Protein contacts are defined by the CaCa distances between residues, and a cutoff
#    of 10 A is used."
DEFAULT_CUTOFF = 10.0

# "giving 1000 appropriately distributed decoys for each contact"
DEFAULT_N_DECOYS = 1000

# Side-chain relaxation protocol. The paper says "short Monte-Carlo relaxation" without
# naming a mover, and the readings differ in cost and outcome (measured on a shuffled
# ubiquitin decoy; 1000-decoy cost in brackets):
#
#   "min"    chi-only gradient minimisation   0.35 s -> 413 REU   [~6 min]
#   "relax"  FastRelax, backbone frozen       3.63 s -> 429 REU   [~1 h]
#
# "relax" is the more literal reading -- FastRelax is Monte-Carlo based, minimisation is
# not. Since the relaxation protocol sets the decoy spread, and that spread is the
# DENOMINATOR of Eq. 1, this choice moves every frustration index. Settled empirically
# in validation, not by preference.
RELAX_PROTOCOLS = ("relax", "min", "none")

# Classification thresholds, taken from frustratometeR so output is comparable.
# NOTE: calibrated on the AWSEM coarse-grained energy function, NOT atomistic REF2015.
# Almost certainly the wrong cut points here; recalibrating is part of validation.
MINIMALLY_FRUSTRATED = 0.78
HIGHLY_FRUSTRATED = -1.0

# Weight on Eq. 2's many-body background term in
#
#     E_ij = e_ij + w * 1/2 (R_i + R_j)
#
# w = 1 is Eq. 2 as the paper writes it (its literal form is 1/2(R_i + R_j), which
# differs from w=1 only by the e_ij term -- measured to be indistinguishable, since
# e_ij ~ -0.7 REU against 1/2(R_i+R_j) ~ -10 REU).
#
# DEFAULT IS 0, a deliberate departure from the paper. At w=1 the index is degenerate:
# validated against frustratometeR on 1UBQ it reports ZERO frustrated contacts where
# frustratometeR reports 30, and 86% of its per-contact signal is reducible to a
# residue-level quantity. At w=0 the pair-specificity (0.44), spread (0.92) and glycine
# behaviour all match frustratometeR closely. See docs/method.md for the sweep.
#
# This does NOT fix the correlation with frustratometeR, which sits at rho ~ 0.15
# regardless of w. It fixes the character of the index, not its agreement.
DEFAULT_BACKGROUND_WEIGHT = 0.0
