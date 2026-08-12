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
