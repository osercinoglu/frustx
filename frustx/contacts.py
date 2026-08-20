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

from frustx.config import CONTACT_ATOMS, DEFAULT_CONTACT_ATOM, DEFAULT_CUTOFF


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
    # Checked BEFORE parsing, and against CONTACT_ATOMS rather than a literal tuple.
    # It used to run after the loop had already indexed res[atom], so an unknown atom
    # name that happened to exist in the structure was silently collected first and only
    # rejected at the end -- and the literal drifted from config.CONTACT_ATOMS for free.
    if atom not in CONTACT_ATOMS:
        raise ValueError(f"atom must be one of {CONTACT_ATOMS}, got {atom!r}")

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

    if not residues:
        raise ValueError(f"No standard amino-acid residues with CA found in {pdb_path}")

    return residues, np.asarray(coords, dtype=np.float64)


# Rows per block in the chunked distance kernels. 256 x n x 3 float64 is ~6 MB at
# n=1000, which keeps the transient off the heap profile without making the Python-level
# loop long enough to matter.
_BLOCK = 256


def pairwise_distances(coords, block=_BLOCK):
    """Euclidean distance matrix, (n, n), computed a block of rows at a time.

    Bit-identical to the full broadcast `sqrt(((a[:,None]-a[None])**2).sum(-1))` -- every
    entry is an independent reduction over three terms, so splitting the rows changes
    nothing about the arithmetic. Verified by test, because "obviously identical" is how
    a last-bit difference gets into a repo that holds serial and parallel runs to
    bit-equality.

    Chunked because the naive one-liner costs 7x the size of its own result: `delta` is
    24n^2 bytes and `delta**2` another 24n^2, against 8n^2 for the `dist` that survives.
    Measured peak was 56 MB at n=1000 and 224 MB at n=2000 -- the comment this replaces
    claimed 8 MB and 200 MB-at-n=5000, i.e. it costed the surviving array and ignored the
    two temporaries. True cost at n=5000 was ~1.4 GB.
    """
    coords = np.asarray(coords, dtype=np.float64)
    n = len(coords)
    out = np.empty((n, n), dtype=np.float64)
    for start in range(0, n, block):
        stop = min(start + block, n)
        delta = coords[start:stop, None, :] - coords[None, :, :]
        out[start:stop] = np.sqrt((delta ** 2).sum(axis=-1))
    return out


def min_heavy_distances(coords, offsets, block=_BLOCK):
    """Minimum heavy-atom distance between every pair of residues, (n, n).

    This is the contact rule a ligand needs. A ligand has no CA and no CB, so the
    CA-CA criterion does not merely give it a bad answer -- it gives it NO answer, and
    `np.nan <= cutoff` is False, which is why contact_pairs rejects non-finite
    coordinates rather than letting a ligand silently touch nothing.

    `coords` is (A, 3), ALL heavy atoms of all residues concatenated; `offsets` is
    (n + 1,) giving each residue's slice, so residues with different atom counts pack
    without padding.

    Chunked over residues because this is the one quantity here that is NOT n^2: with
    ~8 heavy atoms per residue, A = 8n, and a full A x A matrix is 512 MB at n=1000 and
    ~13 GB at n=5000. A naive broadcast is not viable at complex scale even though every
    structure this repo runs today is 76-171 residues.

    The reduction is two segment-minima: atoms -> residues along each axis in turn.
    """
    coords = np.asarray(coords, dtype=np.float64)
    offsets = np.asarray(offsets, dtype=np.intp)
    n = len(offsets) - 1
    if offsets[0] != 0 or offsets[-1] != len(coords):
        raise ValueError(
            f"offsets must span the coordinate array: got {offsets[0]}..{offsets[-1]} "
            f"for {len(coords)} atoms"
        )
    if np.any(np.diff(offsets) < 1):
        # reduceat does NOT return the identity for an empty segment -- it returns the
        # element at that index, so an atomless residue would silently report the
        # distance of whichever atom happened to be next. Refuse instead.
        bad = np.nonzero(np.diff(offsets) < 1)[0]
        raise ValueError(f"residues {bad.tolist()[:5]} have no heavy atoms")

    out = np.empty((n, n), dtype=np.float64)
    starts = offsets[:-1]
    for r0 in range(0, n, block):
        r1 = min(r0 + block, n)
        a0, a1 = offsets[r0], offsets[r1]
        d = np.sqrt((((coords[a0:a1, None, :] - coords[None, :, :]) ** 2).sum(axis=-1)))
        # min over the atoms of each column-residue, then over the atoms of each
        # row-residue in this block
        d = np.minimum.reduceat(d, starts, axis=1)
        out[r0:r1] = np.minimum.reduceat(d, starts[r0:r1] - a0, axis=0)
    return out


class ContactGeometry:
    """Distances computed once per definition; selections are then cheap.

    The motivating waste: scripts/contact_definition.py sweeps cutoffs 10.0/9.5/9.0 over
    seven structures and calls contact_pairs each time, so it pays for 21 full distance
    passes where 7 would do. Nothing about a cutoff changes the distances.

    Definitions are columns rather than an argument to one rule because they are not
    interchangeable: "CA" and "CB" are per-residue representative points, "heavy" is a
    minimum over atom pairs, and a mixed protein/ligand structure needs BOTH at once --
    a CA-CA cutoff between two protein residues and a heavy-atom minimum wherever a
    ligand is involved.

    contact_pairs is deliberately left alone. It takes coordinates and never an atom
    name (tests/test_contact_atom.py pins that), which is the right shape for the
    single-definition case and is what every caller in scripts/ uses today.
    """

    def __init__(self, residues, columns):
        """`columns` maps a definition name to either an (n, 3) coordinate array or a
        precomputed (n, n) distance matrix."""
        self.residues = residues
        n = len(residues)
        self._d = {}
        for name, value in columns.items():
            value = np.asarray(value, dtype=np.float64)
            if value.shape == (n, 3):
                value = pairwise_distances(value)
            elif value.shape != (n, n):
                raise ValueError(
                    f"column {name!r} has shape {value.shape}, expected ({n}, 3) "
                    f"coordinates or ({n}, {n}) distances"
                )
            self._d[name] = value

    @property
    def definitions(self):
        return tuple(sorted(self._d))

    def distances(self, definition):
        if definition not in self._d:
            raise KeyError(
                f"no column {definition!r}; have {self.definitions}"
            )
        return self._d[definition]

    def pairs(self, definition=DEFAULT_CONTACT_ATOM, cutoff=DEFAULT_CUTOFF,
              min_seq_sep=1):
        """The contact set under one definition. Same semantics as contact_pairs."""
        return select_pairs(self.residues, self.distances(definition),
                            cutoff=cutoff, min_seq_sep=min_seq_sep)


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

    return select_pairs(residues, pairwise_distances(ca_coords),
                        cutoff=cutoff, min_seq_sep=min_seq_sep)


def select_pairs(residues, dist, cutoff=DEFAULT_CUTOFF, min_seq_sep=1):
    """Contact pairs from an ALREADY-COMPUTED distance matrix.

    Split out of contact_pairs so that a distance matrix computed once can be selected
    over many times -- different cutoffs, or a different rule per residue kind -- without
    recomputing the geometry. contact_pairs remains the coordinate-taking entry point
    everything in scripts/ uses.

    `dist` may legitimately contain inf (that is how a per-kind rule excludes a pair it
    does not apply to), but NOT NaN: `np.nan <= cutoff` is False, so a NaN would drop the
    pair silently rather than raise. The caller-facing guard against that lives in
    contact_pairs, where the coordinates still exist to be named.
    """
    if min_seq_sep < 1:
        raise ValueError("min_seq_sep must be >= 1 (a residue cannot contact itself)")
    n = len(residues)
    if dist.shape != (n, n):
        raise ValueError(f"dist is {dist.shape}, expected ({n}, {n})")
    if np.isnan(dist).any():
        raise ValueError(
            "distance matrix contains NaN, which would silently drop pairs rather than "
            "raise; use inf to mark a pair the rule does not apply to"
        )

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
