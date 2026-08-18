"""Command-line entry point.

    frustx structure.pdb -o results/

Uses stdlib argparse rather than typer/click: this is a handful of flags and the
dependency would not pay for itself.

Every run writes a run.json alongside the results recording the input, all parameters
and the package version. Frustration indices are only comparable between runs made with
the same decoy count and relaxation protocol, so a results directory that cannot say how
it was produced is not much use.
"""

import argparse
import json
import sys
import time
from pathlib import Path

from frustx import __version__
# Only frustx.config at module scope: it is import-free by design, so `frustx --help`
# and `--version` stay instant. Everything that touches PyRosetta is imported in main().
from frustx.config import (DEFAULT_BACKGROUND_WEIGHT, DEFAULT_CUTOFF,
                          DEFAULT_N_DECOYS, DEFAULT_READOUT, READOUT_SCOPES,
                          CONTACT_ATOMS, DEFAULT_CONTACT_ATOM, RELAX_PROTOCOLS)


def build_parser():
    p = argparse.ArgumentParser(
        prog="frustx",
        description="Local frustration in protein structures at atomistic resolution "
                    "(Chen et al., Nat Commun 11:5944, 2020).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("structure", type=Path, help="input PDB file")
    p.add_argument("-o", "--out", type=Path, required=True,
                   help="output directory (created if absent)")
    p.add_argument("-n", "--decoys", type=int, default=DEFAULT_N_DECOYS,
                   help="number of shuffled decoys")
    p.add_argument("--seed", type=int, default=0,
                   help="base RNG seed for sequence shuffling; decoy k uses seed+k")
    p.add_argument("--protocol", choices=RELAX_PROTOCOLS, default="relax",
                   help="side-chain relaxation: 'relax' is FastRelax (the paper's "
                        "'short Monte-Carlo relaxation'), 'min' is ~10x faster "
                        "gradient minimisation, 'none' skips it")
    p.add_argument("--repeats", type=int, default=1,
                   help="FastRelax repeats; ignored unless --protocol relax")
    p.add_argument("--cutoff", type=float, default=DEFAULT_CUTOFF,
                   help="Ca-Ca contact cutoff in Angstrom")
    p.add_argument("--min-seq-sep", type=int, default=1,
                   help="minimum |i-j| within a chain for a pair to count as a contact; "
                        "1 keeps sequential neighbours, which is the paper's literal reading")
    p.add_argument("--background-weight", type=float, default=DEFAULT_BACKGROUND_WEIGHT,
                   help="weight w on the many-body background in "
                        "E_ij = e_ij + w/2 (R_i + R_j). 0 (default) uses the direct pair "
                        "energy alone; 1 reproduces the paper's Eq. 2, which was measured "
                        "to give a degenerate index -- see docs/method.md")
    p.add_argument("--contact-atom", choices=CONTACT_ATOMS, default=DEFAULT_CONTACT_ATOM,
                   help="atom representing a residue in the contact map. 'CA' (default) "
                        "is the paper's definition; 'CB' (Gly falls back to CA) admits "
                        "far fewer pairs with no atomistic interaction -- see config.py")
    p.add_argument("--readout", choices=READOUT_SCOPES, default=DEFAULT_READOUT,
                   help="energy scope Eq. 1 is applied to: 'pair' (default, the bare "
                        "contact energy) or 'neighbourhood' (that contact plus every "
                        "other one touching i or j, matching what frustratometeR sums). "
                        "Diagnostic -- see docs/method.md before using it")
    p.add_argument("--packing-frustration", action="store_true",
                   help="keep the repulsive fa_rep term when measuring e_ij. This is the "
                        "paper's separate 'packing frustration', which diagnoses structure "
                        "quality rather than functional frustration -- NOT the default index")
    p.add_argument("-j", "--jobs", type=int, default=1,
                   help="decoy-building processes (default 1, serial). Decoys are "
                        "independent, so this scales nearly linearly up to the physical "
                        "core count; measured ~3.6x at -j 4 on 4 physical cores")
    p.add_argument("--packing-seed", type=int, default=None,
                   help="pin decoy k's packing to PACKING_SEED+k, making the ensemble a "
                        "pure function of the seed and identical at any --jobs. Off by "
                        "default: an unseeded packer is what makes a re-run an "
                        "independent sample. Mainly for verification")
    p.add_argument("--jran-base", type=int, default=None,
                   help="base for the per-worker Rosetta RNG streams. Default draws "
                        "fresh entropy per invocation, so re-running with --jobs grows "
                        "the ensemble as it does serially. Set it to reproduce a "
                        "specific parallel run; the resolved value is in run.json")
    p.add_argument("--quiet", action="store_true", help="suppress progress output")
    p.add_argument("--version", action="version", version=f"frustx {__version__}")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    if not args.structure.is_file():
        print(f"error: no such file: {args.structure}", file=sys.stderr)
        return 2
    if args.decoys < 2:
        print("error: --decoys must be at least 2 to estimate a standard deviation",
              file=sys.stderr)
        return 2

    # Imported here, not at module scope: PyRosetta takes ~2.5 s to load, and
    # `frustx --help` must not pay for it. Everything below this line is heavy.
    import pyrosetta
    from frustx.energies import init_rosetta, make_score_function
    from frustx.frustration import compute_frustration
    from frustx.output import contact_table, residue_table, write_bfactor_pdb

    args.out.mkdir(parents=True, exist_ok=True)

    init_rosetta()
    pose = pyrosetta.pose_from_pdb(str(args.structure))

    if not args.quiet:
        print(f"{args.structure.name}: {pose.total_residue()} residues, "
              f"{args.decoys} decoys, protocol={args.protocol}", file=sys.stderr)

    def progress(done, total):
        if not args.quiet:
            # \r keeps this to one line; the final newline is printed after the loop.
            print(f"\r  decoy {done}/{total}", end="", file=sys.stderr, flush=True)

    started = time.time()
    result = compute_frustration(
        pose,
        # Structures are BUILT with fa_rep (it is what stops rotamers overlapping) and
        # MEASURED without it, unless the user asked for packing frustration.
        sf_pack=make_score_function(remove_fa_rep=False),
        sf_measure=make_score_function(remove_fa_rep=not args.packing_frustration),
        n_decoys=args.decoys,
        seed=args.seed,
        protocol=args.protocol,
        repeats=args.repeats,
        cutoff=args.cutoff,
        min_seq_sep=args.min_seq_sep,
        background_weight=args.background_weight,
        readout=args.readout,
        contact_atom=args.contact_atom,
        n_jobs=args.jobs,
        packing_seed=args.packing_seed,
        jran_base=args.jran_base,
        progress=progress,
    )
    elapsed = time.time() - started
    if not args.quiet:
        print(file=sys.stderr)

    contacts = contact_table(result)
    residues = residue_table(result)
    contacts.to_csv(args.out / "contacts.csv", index=False)
    residues.to_csv(args.out / "residues.csv", index=False)
    write_bfactor_pdb(pose.clone(), result, args.out / "frustration.pdb")

    # Provenance: without this a results directory cannot say how it was produced.
    (args.out / "run.json").write_text(json.dumps({
        "frustx_version": __version__,
        "structure": str(args.structure.resolve()),
        "n_residues": pose.total_residue(),
        "n_contacts": int(len(result.contacts)),
        "n_decoys": args.decoys,
        "seed": args.seed,
        "protocol": args.protocol,
        "repeats": args.repeats,
        "cutoff": args.cutoff,
        "min_seq_sep": args.min_seq_sep,
        "background_weight": args.background_weight,
        "readout": args.readout,
        "contact_atom": args.contact_atom,
        "n_jobs": args.jobs,
        "packing_seed": args.packing_seed,
        "jran_base": result.jran_base,
        "packing_frustration": args.packing_frustration,
        "elapsed_seconds": round(elapsed, 1),
    }, indent=2) + "\n")

    if not args.quiet:
        counts = contacts["frustration_class"].value_counts().to_dict()
        print(f"{len(contacts)} contacts in {elapsed:.0f}s -> {args.out}", file=sys.stderr)
        print(f"  {counts}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
