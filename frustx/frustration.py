"""Eq. 1: the local frustration index.

    F_ij = ( <E_decoy> - E_native ) / sigma(E_decoy)

SIGN CONVENTION -- this is the negative of the paper's Eq. 1, deliberately.

The paper writes F = (E0 - <E_U>)/sigma, which makes a *minimally* frustrated contact
come out negative: the native contact energy sits below the decoy mean. Every existing
frustration tool reports the opposite sign, with frustratometeR classifying minimally
frustrated as F > 0.78 and highly frustrated as F < -1. Measured on ubiquitin, the
paper's literal formula gives us a range of -4.05 to +0.55 -- flipped relative to the
convention, as expected. We negate here so that in FrustX output, as everywhere else in
the field:

    high positive F  ->  minimally frustrated (native much better than random)
    near zero        ->  neutral
    negative F       ->  highly frustrated (native no better, or worse, than random)

The contact energy E_ij comes from Eq. 2, which simplifies exactly to 1/2 (R_i + R_j).
See docs/method.md for the derivation and its numerical verification.

KNOWN ARTIFACT -- glycine positions.

Whole-sequence shuffling drops bulky residues onto positions where only glycine fits
sterically. Decoy energies there balloon, and the native glycine looks superb by
comparison, producing a strongly "minimally frustrated" reading. On ubiquitin the
C-terminal LRGG tail scores among the most minimally frustrated contacts
(GLY75-GLY76 at F=+3.66 in our sign) even though that tail is flexible, solvent-exposed
and the functional conjugation site -- it should read as frustrated.

This is a composition-and-sterics effect, not folding frustration. It is inherent to the
paper's shuffling protocol rather than specific to this implementation, and the paper does
not address it. Treat glycine-rich and tail regions with suspicion until the
frustratometeR comparison shows whether the reference implementation behaves the same way.
"""

from dataclasses import dataclass

import numpy as np

from frustx.contacts import DEFAULT_CUTOFF, Residue, contact_pairs
from frustx.decoys import make_decoy, native_reference
from frustx.energies import pair_energy_matrix

from frustx.config import (  # constants live there so the CLI can read
    DEFAULT_N_DECOYS,          # them without importing PyRosetta
    HIGHLY_FRUSTRATED,
    MINIMALLY_FRUSTRATED,
)


@dataclass
class FrustrationResult:
    """Everything needed to interpret or re-analyse a run, without recomputing."""

    residues: list          # list[Residue], 0-based, aligned to matrix indices
    contacts: np.ndarray    # (n_contacts, 2) int, i < j
    native_energy: np.ndarray   # (n, n) E_ij of the native
    decoy_mean: np.ndarray      # (n, n) <E_ij> over decoys
    decoy_std: np.ndarray       # (n, n) sigma(E_ij)
    index: np.ndarray           # (n, n) F_ij, conventional sign
    n_decoys: int


def residues_from_pose(pose):
    """Build Residue records from a pose, preserving PDB chain/number/icode.

    Deliberately reads the pose rather than re-parsing the PDB with Biopython: the
    energies come from the pose, so taking residue identity from anywhere else risks
    the two numbering schemes drifting apart silently.
    """
    info = pose.pdb_info()
    residues = []
    for i in range(1, pose.total_residue() + 1):
        if info is not None:
            chain, resseq, icode = info.chain(i), info.number(i), info.icode(i)
        else:
            # Poses built from sequence have no PDB info.
            chain, resseq, icode = "A", i, " "
        residues.append(
            Residue(
                index=i - 1,
                chain=chain,
                resseq=resseq,
                icode=icode,
                resname=pose.residue(i).name3(),
            )
        )
    return residues


def ca_coords_from_pose(pose):
    """CA coordinates, shape (n_residues, 3), aligned to pose numbering."""
    return np.array(
        [np.array(pose.residue(i).xyz("CA")) for i in range(1, pose.total_residue() + 1)]
    )


def contact_energy_matrix(pose, sf_measure):
    """Eq. 2 contact energies for a pose: E_ij = 1/2 (R_i + R_j).

    R_i is residue i's total pairwise interaction energy. The direct e_ij term cancels
    out of Eq. 2 exactly -- derivation and numerical check in docs/method.md.
    """
    R = pair_energy_matrix(pose, sf_measure).sum(axis=1)
    E = 0.5 * (R[:, None] + R[None, :])
    np.fill_diagonal(E, 0.0)
    return E


def compute_frustration(
    native_pose,
    sf_pack,
    sf_measure,
    n_decoys=DEFAULT_N_DECOYS,
    seed=0,
    protocol="relax",
    repeats=1,
    cutoff=DEFAULT_CUTOFF,
    min_seq_sep=1,
    progress=None,
):
    """Run the full protocol and return a FrustrationResult.

    Parameters
    ----------
    sf_pack : score function WITH fa_rep -- builds decoy structures.
    sf_measure : score function WITHOUT fa_rep -- measures e_ij. See energies.py.
    protocol : passed to relax_sidechains. Defaults to "relax" (FastRelax), the literal
        reading of the paper's "short Monte-Carlo relaxation". "min" is ~10x faster.
    progress : optional callable(done, total) invoked after each decoy. A 1000-decoy run
        takes roughly 12 min with protocol="min" and around an hour with "relax", so a
        caller usually wants some feedback.

    The native reference goes through the same repack-and-relax as the decoys, which the
    paper requires -- see decoys.native_reference.
    """
    if n_decoys < 2:
        raise ValueError("n_decoys must be at least 2 to estimate a standard deviation")

    residues = residues_from_pose(native_pose)
    contacts = contact_pairs(
        residues, ca_coords_from_pose(native_pose), cutoff=cutoff, min_seq_sep=min_seq_sep
    )

    native = native_reference(native_pose, sf_pack, protocol=protocol, repeats=repeats)
    E0 = contact_energy_matrix(native, sf_measure)

    # Accumulate running sums rather than holding 1000 (n x n) matrices: at n=500 that
    # would be 2 GB.
    n = E0.shape[0]
    total = np.zeros((n, n))
    total_sq = np.zeros((n, n))
    for k in range(n_decoys):
        decoy = make_decoy(
            native_pose, sf_pack, seed=seed + k, protocol=protocol, repeats=repeats
        )
        E = contact_energy_matrix(decoy, sf_measure)
        total += E
        total_sq += E * E
        if progress is not None:
            progress(k + 1, n_decoys)

    mean = total / n_decoys
    # Population variance, matching the paper's sigma over the N decoys. Clipped at zero
    # because catastrophic cancellation can push it slightly negative.
    var = np.maximum(total_sq / n_decoys - mean * mean, 0.0)
    std = np.sqrt(var)

    # Sign flipped relative to the paper -- see module docstring.
    with np.errstate(divide="ignore", invalid="ignore"):
        index = np.where(std > 0, (mean - E0) / std, np.nan)

    return FrustrationResult(
        residues=residues,
        contacts=contacts,
        native_energy=E0,
        decoy_mean=mean,
        decoy_std=std,
        index=index,
        n_decoys=n_decoys,
    )


def classify(values):
    """Label frustration indices using the frustratometeR cut points.

    See the caveat on MINIMALLY_FRUSTRATED: these thresholds are AWSEM-calibrated and
    are placeholders until we recalibrate against atomistic output.
    """
    values = np.asarray(values, dtype=float)
    out = np.full(values.shape, "neutral", dtype=object)
    out[values >= MINIMALLY_FRUSTRATED] = "minimally"
    out[values <= HIGHLY_FRUSTRATED] = "highly"
    out[~np.isfinite(values)] = "undefined"
    return out
