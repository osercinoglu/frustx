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

Memory is not a concern: 1UBQ is 76 residues, so 200 decoys is ~9 MB. RAM does become a
wall for big proteins -- the tensor is 8*N*n^2 bytes, so n=500 at N=1000 is 2 GB.

CHECKPOINTING. This script used to build the whole tensor in RAM and write only after the
final decoy, so a job killed at 499/500 lost everything -- which is exactly what happened
to a fig2 run. It now writes <out>.partial.npz every CHECKPOINT_EVERY decoys and resumes
from it, re-entering the loop at the first unfinished decoy.

Resumption reproduces the SEQUENCE ensemble exactly -- decoy k's shuffle depends only on
seed=k (frustx/decoys.py:130) -- but NOT the energies bit-for-bit. Rosetta is initialised
without -constant_seed (frustx/energies.py:37), so the packer draws from a global RNG that
differs between processes: decoy k built after a resume has the same sequence as decoy k
built before it, and a slightly different packing. A resumed run is therefore statistically
equivalent to an uninterrupted one, not identical to it. That is the correct standard here
(the decoy ensemble is a sample, not a fixed object), but it does mean a resumed run cannot
be used to reproduce an earlier run's exact numbers.

The native reference is stored in the checkpoint and NOT recomputed on resume. It goes
through a repack whose packer draws from Rosetta's global RNG, which is not seeded per
call, so recomputing it would silently give a different E0 for the second half of the run.

PARALLEL. Decoys are independent, so this forks a worker pool. Two things change from
the serial version and both matter:

  * Completion is OUT OF ORDER, so the checkpoint cannot be a "first n_done decoys are
    finished" high-water mark any more -- that would be meaningless. It is a boolean
    `done` mask of length n_decoys, and a resume submits exactly the unfinished indices.
  * The tensor is allocated with np.full(nan) rather than np.empty. With a prefix count
    and np.empty, a checkpoint claiming n_done=35 while decoys 25-29 were still in flight
    would have shipped UNINITIALISED MEMORY into the final .npz, and nothing downstream
    would have noticed.

Usage:
    .venv/bin/python scripts/dump_decoy_samples.py <pdb> <out.npz> [n_decoys] [protocol] [n_jobs]

    Re-running the same command after an interruption resumes; delete <out>.partial.npz
    to force a clean start.
"""
import os
import sys

import numpy as np
import pyrosetta

from frustx.decoys import make_decoy, native_reference
from frustx.energies import init_rosetta, make_score_function
from frustx.config import DEFAULT_BACKGROUND_WEIGHT
from frustx.frustration import contact_energy_matrix, residues_from_pose


# Every 25 decoys is ~50 s of work at the measured 'min' cost on a 169-residue protein --
# small enough that a teardown costs little, large enough that rewriting the partial
# tensor stays a rounding error on the total.
CHECKPOINT_EVERY = 25


def _load_checkpoint(path, n_decoys, n):
    """Return (decoys, native, done_mask) from a partial run, or (None, None, None).

    A checkpoint whose shape does not match what was asked for is ignored rather than
    trusted -- resuming a 500-decoy run into a 1000-decoy request would silently mix two
    different requests.
    """
    if not os.path.exists(path):
        return None, None, None
    try:
        d = np.load(path, allow_pickle=True)
        done = d["done"].astype(bool)
        if d["decoys"].shape != (n_decoys, n, n) or done.shape != (n_decoys,):
            print(f"ignoring {path}: shape mismatch", flush=True)
            return None, None, None
        print(f"resuming from {path}: {done.sum()}/{n_decoys} decoys already done",
              flush=True)
        return d["decoys"], d["native"], done
    except Exception as exc:                       # a truncated write is not fatal
        print(f"ignoring unreadable checkpoint {path}: {exc}", flush=True)
        return None, None, None


def main(pdb, out, n_decoys, protocol, n_jobs=1):
    init_rosetta()
    pose = pyrosetta.pose_from_pdb(pdb)

    # Same two score functions the production path uses: fa_rep ON for packing so
    # decoys are physically buildable, OFF for measurement so a single clash cannot
    # dominate the contact energy.
    sf_pack = make_score_function(remove_fa_rep=False)
    sf_measure = make_score_function(remove_fa_rep=True)

    n = pose.total_residue()
    partial = str(out) + ".partial.npz"
    decoys, native, done = _load_checkpoint(partial, n_decoys, n)

    if decoys is None:
        # The native goes through the same repack as the decoys -- comparing an
        # unrelaxed native against relaxed decoys would bias every index positive. It is
        # built ONCE, in the parent, and stored in the checkpoint: a worker rebuilding it
        # under a different RNG stream would give a different E0 for part of the run.
        nat_pose = native_reference(pose, sf_pack, protocol=protocol)
        native = contact_energy_matrix(nat_pose, sf_measure)
        # NaN, not np.empty: an unfinished slot must be obviously unfinished. See the
        # module docstring on why np.empty plus a prefix count was actively dangerous.
        decoys = np.full((n_decoys, n, n), np.nan, dtype=np.float64)
        done = np.zeros(n_decoys, dtype=bool)

    res = residues_from_pose(pose)
    labels = dict(
        resseq=np.array([r.resseq for r in res]),
        chain=np.array([r.chain for r in res]),
        resname=np.array([r.resname for r in res]),
    )

    todo = np.flatnonzero(~done).tolist()
    n_jobs = max(1, min(n_jobs, len(todo))) if todo else 1

    def save_checkpoint():
        # Write to a temp path and rename: a kill DURING the checkpoint write would
        # otherwise leave a truncated file where a good one used to be.
        tmp = partial + ".tmp.npz"
        np.savez_compressed(tmp, decoys=decoys, native=native, done=done, **labels)
        os.replace(tmp, partial)

    def record(k, E):
        # Indexed by the k carried in the RESULT, never by arrival order -- under
        # out-of-order completion those are different, and using arrival order would
        # permute the tensor rows. Nothing downstream would complain: mean, sigma and
        # index are all order-invariant, but scripts/packer_noise.py and the
        # Rao-Blackwell analysis join row k to seed k and would be silently wrong.
        decoys[k] = E
        done[k] = True

    print(f"{len(todo)} decoys to build on {n_jobs} worker(s)", flush=True)
    if n_jobs == 1:
        for count, k in enumerate(todo, 1):
            record(k, contact_energy_matrix(
                make_decoy(pose, sf_pack, seed=k, protocol=protocol), sf_measure))
            if count % 10 == 0:
                print(f"{int(done.sum())}/{n_decoys}", flush=True)
            if count % CHECKPOINT_EVERY == 0 and count < len(todo):
                save_checkpoint()
                print(f"  checkpoint at {int(done.sum())}", flush=True)
    else:
        import multiprocessing as mp

        from frustx.frustration import _WORKER, _worker_init, _worker_decoy, contact_mask
        ctx = mp.get_context("fork")
        # readout="pair" with an all-True mask: this script wants the raw contact energy
        # matrix, the same quantity the serial version stored.
        _WORKER.update(pose=pose, sf_pack=sf_pack, sf_measure=sf_measure,
                       mask=None, seed=0, protocol=protocol, repeats=1,
                       background_weight=DEFAULT_BACKGROUND_WEIGHT, readout="pair",
                       packing_seed=None)
        counter = ctx.Value("i", 0)
        try:
            with ctx.Pool(n_jobs, initializer=_worker_init,
                          initargs=(1, counter)) as pool:
                for count, (k, E) in enumerate(
                        pool.imap_unordered(_worker_decoy, todo, chunksize=1), 1):
                    record(k, E)
                    if count % 10 == 0:
                        print(f"{int(done.sum())}/{n_decoys}", flush=True)
                    if count % CHECKPOINT_EVERY == 0 and count < len(todo):
                        save_checkpoint()
                        print(f"  checkpoint at {int(done.sum())}", flush=True)
        finally:
            _WORKER.clear()

    if not done.all():
        raise RuntimeError(f"only {int(done.sum())}/{n_decoys} decoys completed")
    assert not np.isnan(decoys).any(), "a decoy slot was never written"

    np.savez_compressed(out, decoys=decoys, native=native, **labels)
    if os.path.exists(partial):
        os.remove(partial)          # only once the real output is safely on disk
    print(f"wrote {out}  decoys={decoys.shape}")


if __name__ == "__main__":
    # Default protocol is "min" so this matches results/validation/frustx_min500,
    # the run whose index distribution motivated the question.
    main(sys.argv[1], sys.argv[2],
         int(sys.argv[3]) if len(sys.argv) > 3 else 200,
         sys.argv[4] if len(sys.argv) > 4 else "min",
         int(sys.argv[5]) if len(sys.argv) > 5 else 1)
