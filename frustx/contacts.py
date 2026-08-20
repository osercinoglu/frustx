"""Contact-map construction.

The Chen et al. (2020) protocol defines a protein contact purely geometrically:

    "Protein contacts are defined by the CalphaCalpha distances between
     residues, and a cutoff of 10 A is used."

Everything downstream in FrustX is computed *per contact* -- the decoy energy
ensembles, the Eq. 2 many-body contact energies, and the Eq. 1 frustration
indices -- so this module is what fixes the meaning of "a contact" for the whole
package.  Keeping it separate (and free of any PyRosetta dependency) means the
contact definition can be unit-tested on hand-built coordinates.
"""

from dataclasses import dataclass

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

from frustx.config import DEFAULT_CONTACT_ATOM, DEFAULT_CUTOFF


@dataclass(frozen=True)
class Residue:
    """One residue that survived parsing and carries a CA atom.

    `index` is the 0-based position in the parsed list and is the identifier
    used everywhere else in FrustX -- it indexes rows of the coordinate array
    and, later, rows/columns of the Rosetta pair-energy matrix.  It is
    deliberately NOT the PDB residue number, which can restart per chain, skip
    values, or carry insertion codes.
    """

    index: int
    chain: str
    resseq: int
    icode: str
    resname: str

    @property
    def label(self) -> str:
        """Human-readable id, e.g. 'A:LYS745' (or 'A:LYS745A' with an icode)."""
        return f"{self.chain}:{self.resname}{self.resseq}{self.icode.strip()}"


def load_ca(pdb_path, model_id=0, atom=DEFAULT_CONTACT_ATOM):
    """Parse a PDB file and return (residues, coords).

    `atom` selects which atom represents a residue for contact purposes:

        "CA"  the paper's definition (Chen et al. specify Calpha-Calpha).
        "CB"  the side-chain direction, falling back to CA for glycine (which has
              no CB) and for residues whose CB is unresolved. This is
              frustratometeR's convention, and it admits far fewer energetically
              dead pairs -- see DEFAULT_CONTACT_ATOM in config.py for the numbers.

    Returns
    -------
    residues : list[Residue]
    ca_coords : np.ndarray, shape (n_residues, 3), float64

    Only standard amino acids with a CA atom are kept; waters, ligands, ions and
    non-standard residues are dropped.  Only one model is read (NMR ensembles
    would otherwise duplicate every residue).  Where alternate locations exist we
    take whichever altloc Biopython considers the primary one, which is its
    default behaviour.
    """
    # QUIET=1 suppresses Biopython's warnings about discontinuous chains, which
    # are common and harmless in crystal structures with unresolved loops.
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", str(pdb_path))
    model = structure[model_id]

    residues, coords = [], []
    for chain in model:
        for res in chain:
            # is_aa(standard=True) rejects HETATM residues, waters and modified
            # amino acids.  We additionally require a CA: residues with only
            # backbone N/C resolved would otherwise have no coordinate.
            if not is_aa(res, standard=True):
                continue
            if "CA" not in res:
                continue
            hetflag, resseq, icode = res.get_id()
            residues.append(
                Residue(
                    index=len(residues),
                    chain=chain.get_id(),
                    resseq=resseq,
                    icode=icode,
                    resname=res.get_resname(),
                )
            )
            # CB falls back to CA for glycine, which has none, and for residues
            # whose CB is unresolved in the density. Silently, because that is the
            # standard convention -- frustratometeR does the same.
            coords.append(res[atom].get_coord() if atom in res
                          else res["CA"].get_coord())

    if atom not in ("CA", "CB"):
        raise ValueError(f"atom must be 'CA' or 'CB', got {atom!r}")
    if not residues:
        raise ValueError(f"No standard amino-acid residues with CA found in {pdb_path}")

    return residues, np.asarray(coords, dtype=np.float64)


def contact_pairs(residues, ca_coords, cutoff=DEFAULT_CUTOFF, min_seq_sep=1):
    """Return all residue pairs in contact, as an (n_contacts, 2) int array.

    Parameters
    ----------
    cutoff : float
        Ca-Ca distance cutoff in Angstrom.  10.0 in the paper.
    min_seq_sep : int
        Minimum |i - j| in *residue index* for a pair to count, applied only
        within a chain.  min_seq_sep=1 keeps every pair including sequential
        neighbours, which is the literal reading of the paper -- its Methods
        state a distance cutoff and nothing else.

        NOTE (open question, see docs/method.md): other frustration tools
        exclude near-sequence neighbours, whose contact is a foregone conclusion
        of the covalent geometry rather than a packing choice.  We expose the
        knob and settle its value empirically in the frustratometeR validation
        (Phase 6) instead of guessing now.

    Pairs are returned with i < j, sorted, so each contact appears exactly once.
    """
    if min_seq_sep < 1:
        raise ValueError("min_seq_sep must be >= 1 (a residue cannot contact itself)")

    n = len(residues)
    if n != len(ca_coords):
        raise ValueError(f"residues ({n}) and ca_coords ({len(ca_coords)}) disagree")

    ca_coords = np.asarray(ca_coords, dtype=np.float64)

    # A non-finite coordinate is the one input this function cannot fail loudly on by
    # itself: `np.nan <= cutoff` is False, with no warning, so the pair is simply
    # ABSENT from the contact set and the run goes on to report a clean answer over an
    # incomplete map.  Nothing downstream can distinguish that from a residue that
    # genuinely contacts nothing, so it has to be caught here or not at all.
    #
    # The way this arrives is a residue with no representative atom -- a ligand, a
    # water, an ion -- given a NaN placeholder by whatever built the coordinates.  That
    # is the natural shape of a ligand-aware coordinate builder, so this guard is what
    # stops "we forgot to give ligands a contact rule" from looking like "the ligand
    # touches nothing".  See docs/method.md, "What the EGFR atomfrust branch already
    # solved", where exactly this cost a run.
    bad = np.nonzero(~np.isfinite(ca_coords).all(axis=1))[0]
    if bad.size:
        shown = ", ".join(residues[k].label for k in bad[:5])
        more = f" (and {bad.size - 5} more)" if bad.size > 5 else ""
        raise ValueError(
            f"non-finite contact coordinate for {bad.size} residue(s): {shown}{more}. "
            f"A NaN coordinate would silently drop every pair involving these residues "
            f"rather than raising, so it is rejected here."
        )

    # Full pairwise distance matrix by broadcasting.  This is O(n^2) memory
    # (~8 MB at n=1000, ~200 MB at n=5000).  Fine for single chains and typical
    # complexes; if we ever hit ribosome-scale input this is the line to swap
    # for a neighbour-list / KD-tree.
    delta = ca_coords[:, None, :] - ca_coords[None, :, :]
    dist = np.sqrt((delta**2).sum(axis=-1))

    # Upper triangle only (k=1 excludes the diagonal), so each pair is unique.
    close = np.triu(dist <= cutoff, k=1)

    if min_seq_sep > 1:
        # Sequence separation is only meaningful within one chain: residues in
        # different chains are never "sequential neighbours" however close their
        # indices happen to be in our flat numbering.
        chain_ids = np.array([r.chain for r in residues])
        same_chain = chain_ids[:, None] == chain_ids[None, :]
        idx = np.arange(n)
        sep = np.abs(idx[:, None] - idx[None, :])
        too_close_in_sequence = same_chain & (sep < min_seq_sep)
        close &= ~too_close_in_sequence

    i, j = np.nonzero(close)
    return np.column_stack([i, j])
