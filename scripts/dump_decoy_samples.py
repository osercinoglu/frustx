"""Dump the full per-decoy contact energy tensor, not just its mean and sd.

`compute_frustration` accumulates running sums and throws each decoy away, which is
right for production but leaves us unable to ask what *shape* the decoy energy
distribution has. That shape is the whole question behind the classification
thresholds: FrustX's index is a genuine Z-score, so a cut at Z = 1 only carries a
probabilistic meaning ("fewer than ~16% of random substitutions do this well") if the
decoy energies are approximately normal. Borrowing AWSEM's 0.78 / -1.0 assumes a
distribution shape we have never checked.

Output is an .npz holding:
    decoys  (n_decoys, n_res, n_res)  E_ij for every decoy
    native  (n_res, n_res)            E_ij for the repacked native reference
    resseq, chain, resname            residue labels, so contacts can be joined back

Memory is not a concern: 1UBQ is 76 residues, so 200 decoys is ~9 MB.

Usage:
    .venv/bin/python scripts/dump_decoy_samples.py <pdb> <out.npz> [n_decoys] [protocol]
"""
import sys

import numpy as np
import pyrosetta

from frustx.decoys import make_decoy, native_reference
from frustx.energies import init_rosetta, make_score_function
from frustx.frustration import contact_energy_matrix, residues_from_pose


def main(pdb, out, n_decoys, protocol):
    init_rosetta()
    pose = pyrosetta.pose_from_pdb(pdb)

    # Same two score functions the production path uses: fa_rep ON for packing so
    # decoys are physically buildable, OFF for measurement so a single clash cannot
    # dominate the contact energy.
    sf_pack = make_score_function(remove_fa_rep=False)
    sf_measure = make_score_function(remove_fa_rep=True)

    # The native goes through the same repack as the decoys -- comparing an
    # unrelaxed native against relaxed decoys would bias every index positive.
    nat_pose = native_reference(pose, sf_pack, protocol=protocol)
    native = contact_energy_matrix(nat_pose, sf_measure)

    n = pose.total_residue()
    decoys = np.empty((n_decoys, n, n), dtype=np.float64)
    for k in range(n_decoys):
        # Seed by index so the run is reproducible and resumable.
        decoys[k] = contact_energy_matrix(make_decoy(pose, sf_pack, seed=k, protocol=protocol), sf_measure)
        if (k + 1) % 10 == 0:
            print(f"{k + 1}/{n_decoys}", flush=True)

    res = residues_from_pose(pose)
    np.savez_compressed(
        out,
        decoys=decoys,
        native=native,
        resseq=np.array([r.resseq for r in res]),
        chain=np.array([r.chain for r in res]),
        resname=np.array([r.resname for r in res]),
    )
    print(f"wrote {out}  decoys={decoys.shape}")


if __name__ == "__main__":
    # Default protocol is "min" so this matches results/validation/frustx_min500,
    # the run whose index distribution motivated the question.
    main(sys.argv[1], sys.argv[2],
         int(sys.argv[3]) if len(sys.argv) > 3 else 200,
         sys.argv[4] if len(sys.argv) > 4 else "min")
