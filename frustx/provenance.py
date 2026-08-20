"""Content keys for deciding when previously computed work is still valid.

A run is a chain of stages, and they do NOT all depend on the same settings. Changing
`--readout` invalidates how decoys are MEASURED but not the decoy structures themselves,
which are the part that costs ~700 ms each. Changing `--protocol` invalidates both. A
single "settings changed" flag cannot express that, and the alternative -- checking a few
fields by hand at each reuse point -- is what this module exists to replace, because the
hand-written version silently omits whatever was added to the CLI last.

The immediate consumer is checkpoint resumption in scripts/dump_decoy_samples.py, where
the pre-existing check was a SHAPE comparison: it caught a changed decoy count or a
changed residue count and nothing else, so resuming a `min` checkpoint under `relax`
mixed two ensembles into one tensor while the native reference stayed pinned to the
first protocol. See "Tier-2 ports" in docs/method.md.

Deliberately NOT a general caching layer. Nothing here decides what to recompute; it
only answers "were these produced under the same settings as those?".
"""

import hashlib
import json
from pathlib import Path

# What each stage actually reads. Taken from tracing compute_frustration, not from the
# CLI's argument list -- the two differ, and it is the trace that is correct.
#
#   structure    parsing the input into a pose
#   contacts     the contact map
#   ensemble     BUILDING decoy structures (the expensive stage)
#   measurement  reading energies off those structures (cheap, but a superset of
#                `ensemble`'s settings -- which is exactly why they are separate)
#
# `weights` and `init_flags` appear in several stages because a score function and a
# residue type set are inputs to all of them, even though neither is a CLI flag today.
STAGE_INPUTS = {
    "structure": ("structure_key", "init_flags"),
    "contacts": ("structure_key", "init_flags", "contact_atom", "cutoff", "min_seq_sep"),
    "ensemble": ("structure_key", "init_flags", "weights", "protocol", "repeats",
                 "seed", "n_decoys", "packing_seed"),
    "measurement": ("structure_key", "init_flags", "weights", "protocol", "repeats",
                    "seed", "n_decoys", "packing_seed",
                    "background_weight", "readout", "packing_frustration"),
}


def file_key(path, n_bytes=1 << 20):
    """A content hash of an input file, NOT its path.

    Two runs against different edits of the same filename are otherwise
    indistinguishable, which is the failure mode a provenance record most needs to
    catch -- a path proves nothing about what was in the file at the time.

    Read in chunks because a structure file is small today but this should not be the
    line that fails on a large complex.
    """
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(n_bytes), b""):
            h.update(chunk)
    return h.hexdigest()


def regeneration_key(settings, stage):
    """A hex digest over exactly the settings `stage` depends on.

    Equal keys mean the two runs agree on everything that stage reads. They do NOT mean
    the bytes will match: by default this project draws fresh RNG entropy per invocation
    (see the jran_base rationale in frustration.py), so an ensemble is reproducible on
    demand, with --packing-seed, rather than by accident. A key answers "are these
    comparable", not "are these identical".

    `readout` is a special case worth stating: at the default readout="pair" the contact
    map never enters the energies -- apply_readout returns E untouched -- so cutoff,
    min_seq_sep and contact_atom do not belong in the measurement key. At
    readout="neighbourhood" they do. Resolving that here rather than in a fixed table is
    the difference between a key that over-invalidates and one that is correct.
    """
    if stage not in STAGE_INPUTS:
        raise ValueError(f"unknown stage {stage!r}, expected one of {tuple(STAGE_INPUTS)}")

    fields = list(STAGE_INPUTS[stage])
    if stage == "measurement" and settings.get("readout") != "pair":
        fields += ["contact_atom", "cutoff", "min_seq_sep"]

    missing = [f for f in fields if f not in settings]
    if missing:
        # Raise rather than treat absent as a default. A key computed over a silently
        # missing field would COMPARE EQUAL to one computed with the field present at
        # its default, which is the one way this module could actively cause the bug it
        # exists to prevent.
        raise KeyError(f"stage {stage!r} needs settings {missing} and they were not given")

    # sort_keys so the digest does not depend on dict insertion order; default=str so a
    # Path or a numpy scalar does not raise here, far from where it was introduced.
    payload = json.dumps({f: settings[f] for f in sorted(fields)},
                         sort_keys=True, default=str)
    return hashlib.sha256(payload.encode()).hexdigest()


def describe_mismatch(old, new, stage):
    """Which settings differ, for an error message a human can act on.

    A key comparison can only say yes or no; being told "the checkpoint does not match"
    without being told that `protocol` went from min to relax is not enough to decide
    whether to delete the file or fix the command line.
    """
    fields = list(STAGE_INPUTS[stage])
    if stage == "measurement" and (old.get("readout") != "pair"
                                   or new.get("readout") != "pair"):
        fields += ["contact_atom", "cutoff", "min_seq_sep"]
    out = []
    for f in sorted(set(fields)):
        a, b = old.get(f, "<absent>"), new.get(f, "<absent>")
        if a != b:
            out.append(f"{f}: {a!r} -> {b!r}")
    return out
