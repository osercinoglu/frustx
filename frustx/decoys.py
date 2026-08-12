"""Decoy generation: the randomised reference ensemble for Eq. 1.

The paper's protocol, verbatim:

    "we randomly shuffle the protein sequence and then repack the resulting sequence
     onto the backbone that is provided without perturbing the backbone coordinates
     within each protein chain to make sure that only the side chains are re-packed.
     A short Monte-Carlo relaxation is then performed to better eliminate many of the
     possible side-chain clashes with the backbone fixed."

and, critically, for the native reference:

    "The contact energies of native sequence (E0) are obtained in a similar fashion by
     omitting the shuffling step."

That last sentence is not a detail.  The native pose must go through the *same*
repack-and-relax as the decoys, otherwise E0 and the decoy ensemble are measured on
structures prepared differently and the Z-score of Eq. 1 compares two incomparable
things.  `native_reference()` exists solely to make that hard to get wrong.

Why shuffling rather than random mutation: a permutation preserves the native amino-acid
composition exactly, which is what satisfies the paper's "sequence space is randomly
sampled according to the native amino acid frequency distribution".

Two score functions are in play and they are not interchangeable:

  * `sf_pack` KEEPS fa_rep.  It builds the structure, and fa_rep is the only thing
    stopping rotamers from occupying the same space.  Packing without it produces
    physically meaningless overlapping side chains.
  * `sf_measure` DROPS fa_rep.  It is used later, in energies.py, to read e_ij.

Measured on ubiquitin: backbone atoms move by exactly 0.0 A through both design and
relaxation, so "the backbone is fixed" is enforced, not merely intended.
"""

import random

import numpy as np

import pyrosetta
from pyrosetta.rosetta.core.chemical import aa_from_oneletter_code
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.protocols.minimization_packing import MinMover, PackRotamersMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.utility import vector1_bool

from frustx.config import RELAX_PROTOCOLS  # see there for measured costs


def _sidechain_only_movemap():
    """MoveMap allowing side-chain torsions only. This is what freezes the backbone."""
    mm = MoveMap()
    mm.set_bb(False)
    mm.set_chi(True)
    return mm


def shuffle_sequence(sequence, rng):
    """Return a random permutation of `sequence`.

    Composition is preserved exactly, which is the point -- see module docstring.
    `rng` is an explicit random.Random so decoy generation is reproducible.
    """
    chars = list(sequence)
    rng.shuffle(chars)
    return "".join(chars)


def apply_sequence(pose, sequence, sf_pack):
    """Design `pose` to `sequence` in place, repacking side chains on a fixed backbone.

    Implemented by restricting each position to exactly one allowed amino acid, which
    turns Rosetta's design machinery into "mutate to this sequence, then pack".
    """
    if len(sequence) != pose.total_residue():
        raise ValueError(
            f"sequence length {len(sequence)} != pose residues {pose.total_residue()}"
        )

    task = TaskFactory.create_packer_task(pose)
    for i, aa in enumerate(sequence, start=1):
        # A 20-long bool vector, one flag per canonical amino acid.
        keep = vector1_bool(20)
        keep[int(aa_from_oneletter_code(aa))] = True
        task.nonconst_residue_task(i).restrict_absent_canonical_aas(keep)

    PackRotamersMover(sf_pack, task).apply(pose)
    return pose


def repack(pose, sf_pack):
    """Repack side chains without changing the sequence, backbone fixed.

    This is the native reference's counterpart to `apply_sequence`: same packing, no
    sequence change.
    """
    task = TaskFactory.create_packer_task(pose)
    task.restrict_to_repacking()
    PackRotamersMover(sf_pack, task).apply(pose)
    return pose


def relax_sidechains(pose, sf_pack, protocol="relax", repeats=1):
    """Relieve side-chain clashes with the backbone held fixed.

    See RELAX_PROTOCOLS for what the choice costs and why it is not settled.
    """
    if protocol not in RELAX_PROTOCOLS:
        raise ValueError(f"protocol must be one of {RELAX_PROTOCOLS}, got {protocol!r}")

    if protocol == "none":
        return pose

    mm = _sidechain_only_movemap()
    if protocol == "min":
        MinMover(mm, sf_pack, "lbfgs_armijo_nonmonotone", 1e-2, True).apply(pose)
    else:
        mover = FastRelax(sf_pack, repeats)
        mover.set_movemap(mm)
        mover.apply(pose)
    return pose


def make_decoy(native_pose, sf_pack, seed, protocol="relax", repeats=1):
    """Build one decoy: shuffle -> design onto fixed backbone -> relax.

    `native_pose` is never modified; a clone is returned.
    """
    rng = random.Random(seed)
    sequence = shuffle_sequence(native_pose.sequence(), rng)

    pose = native_pose.clone()
    apply_sequence(pose, sequence, sf_pack)
    relax_sidechains(pose, sf_pack, protocol=protocol, repeats=repeats)
    return pose


def native_reference(native_pose, sf_pack, protocol="relax", repeats=1):
    """Build the native reference for E0: repack -> relax, no shuffling.

    Deliberately mirrors `make_decoy` step for step. E0 and the decoy energies must come
    off structures prepared identically, or the Z-score is meaningless.
    """
    pose = native_pose.clone()
    repack(pose, sf_pack)
    relax_sidechains(pose, sf_pack, protocol=protocol, repeats=repeats)
    return pose


def backbone_coords(pose):
    """N, CA, C coordinates, shape (3 * n_residues, 3).

    Used by the tests to assert the backbone never moves. Cheap enough to call in
    production as a safety check.
    """
    return np.array(
        [
            np.array(pose.residue(i).xyz(atom))
            for i in range(1, pose.total_residue() + 1)
            for atom in ("N", "CA", "C")
        ]
    )
