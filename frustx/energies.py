"""Pairwise residue interaction energies e_ij from the Rosetta REF2015 force field.

This module answers one question: for a given pose, what is the interaction energy
e_ij between every pair of residues?  Eq. 2 of the paper builds its many-body contact
energy out of those numbers, so getting the extraction exactly right matters more than
anything downstream.

Two non-obvious things are required to make the decomposition complete.  Both were
established empirically -- see "Implementation notes" in docs/method.md for the
measurements.

1. Backbone-backbone hydrogen bonds are NOT written to the energy graph by default.
   On ubiquitin that silently hides 41.96 REU of hbond_sr_bb + hbond_lr_bb -- exactly
   the interactions holding secondary structure together.  Fixed by an option on the
   score function (NOT a command-line flag, which does not exist).

2. rama_prepro is not on the energy graph either; it lives in a long-range container as
   strict (i, i+1) pairs.  It cannot be dismissed as constant under sequence shuffling:
   with the backbone held fixed, single mutations move it by 2.6 to 12.5 REU.

With both handled, one-body + graph edges + long-range containers reproduces the REF2015
total to machine precision.  `tests/test_energies.py` asserts exactly that, so a future
Rosetta version that moves a term elsewhere fails loudly instead of quietly biasing
every frustration index.
"""

import numpy as np

import pyrosetta
from pyrosetta.rosetta.core.scoring import EMapVector, ScoreType, methods

# REF2015: "we employed the REF2015 version of the rosetta energy function, which has a
# set of well-tested weights for each energy term."
DEFAULT_WEIGHTS = "ref2015"

# Rosetta is chatty and aborts on unrecognised residues; neither is useful here.
DEFAULT_INIT_FLAGS = "-mute all -ignore_unrecognized_res"

_INITIALISED = False


def init_rosetta(extra_flags=""):
    """Initialise PyRosetta once per process. Safe to call repeatedly.

    Rosetta's init is global and expensive, and re-initialising mid-run can silently
    reset options, so we guard it.
    """
    global _INITIALISED
    if _INITIALISED:
        return
    pyrosetta.init(f"{DEFAULT_INIT_FLAGS} {extra_flags}".strip(), silent=True)
    _INITIALISED = True


def make_score_function(remove_fa_rep=True, weights=DEFAULT_WEIGHTS):
    """Build the REF2015 score function configured for frustration analysis.

    Parameters
    ----------
    remove_fa_rep : bool
        Zero the repulsive Lennard-Jones weight, as the paper requires:

            "it is both easy and appropriate to simply remove the harsh rapidly
             varying repulsive force term."

        Their reasoning is that residual clashes survive relaxation and inflate the
        decoy variance, which is the denominator of the frustration index -- an
        inflated sigma washes out real signal.  Setting this False instead yields what
        the authors call "packing frustration", a different (also useful) quantity that
        diagnoses structure quality rather than functional frustration.
    """
    sf = pyrosetta.create_score_function(weights)

    # Push backbone-backbone hbonds onto energy-graph edges.  Without this they are
    # accumulated as whole-structure energy and are invisible to a pairwise read.
    # Verified: the weighted total is identical with this on or off -- it changes where
    # energy is stored, not the physics.
    opts = sf.energy_method_options()
    opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    sf.set_energy_method_options(opts)

    if remove_fa_rep:
        sf.set_weight(ScoreType.fa_rep, 0.0)

    return sf


def pair_energy_matrix(pose, sf):
    """Weighted pairwise interaction energies e_ij.

    Returns a symmetric (n_residues, n_residues) float array in REU with a zero
    diagonal, indexed 0-based (Rosetta numbers residues from 1).

    Each entry is the *weighted* sum over all score terms of the interaction between
    that pair -- i.e. e_ij as Eq. 2 uses it.  One-body terms (fa_dun, ref, p_aa_pp,
    omega, fa_intra_*) are deliberately absent: they are properties of a single residue,
    not interactions between two, and Eq. 2 is built purely from pair terms.  A
    consequence worth remembering: sum(e_ij) does NOT equal the pose total score.
    """
    sf(pose)  # scoring is what populates pose.energies(); without it the graph is stale
    energies = pose.energies()
    weights = sf.weights()
    n = pose.total_residue()

    e = np.zeros((n, n), dtype=np.float64)

    # --- short-range pair terms, from the energy graph ----------------------
    # Iterate the edge list directly rather than probing all n^2 pairs with
    # find_energy_edge: only ~1600 edges exist for a 76-residue protein, and this has to
    # run once per decoy.  EnergyEdge.dot(weights) does the weighted sum in C++.
    graph = energies.energy_graph()
    it, end = graph.const_edge_list_begin(), graph.const_edge_list_end()
    while it.__ne__(end):
        edge = it.dereference()
        i = edge.get_first_node_ind() - 1
        j = edge.get_second_node_ind() - 1
        val = edge.dot(weights)
        e[i, j] += val
        e[j, i] += val
        it.pre_increment()

    # --- long-range pair terms, from their containers -----------------------
    # Swept generically rather than hard-coding ramaprepro_lr: dslf_fa13 appears the same
    # way for disulfide-bonded proteins, and our ubiquitin test case has no disulfides,
    # so hard-coding would have passed the test and been wrong on the next protein.
    for lr_type in methods.LongRangeEnergyType.__members__.values():
        container = energies.long_range_container(lr_type)
        if container is None or container.empty():
            continue
        for res in range(1, n + 1):
            it = container.const_upper_neighbor_iterator_begin(res)
            end = container.const_upper_neighbor_iterator_end(res)
            while it.__ne__(end):
                emap = EMapVector()
                it.retrieve_energy(emap)
                i = res - 1
                j = it.upper_neighbor_id() - 1
                val = emap.dot(weights)
                e[i, j] += val
                e[j, i] += val
                it.pre_increment()

    return e


def onebody_energy_vector(pose, sf):
    """Weighted one-body energy per residue, shape (n_residues,), 0-based.

    Not used by Eq. 2 -- these are single-residue properties, not interactions.  It
    exists so the completeness test can verify that one-body + pairwise accounts for the
    entire score, which is the guarantee that no pair energy is being dropped.
    """
    sf(pose)
    energies = pose.energies()
    weights = sf.weights()
    n = pose.total_residue()
    return np.array(
        [energies.onebody_energies(i).dot(weights) for i in range(1, n + 1)],
        dtype=np.float64,
    )
