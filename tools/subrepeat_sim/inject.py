#!/usr/bin/env python3
"""Generate a small synthetic corpus that tests the rescore subrepeat
metric against controlled cases of nested tandem repeats.

Each case starts from a randomly generated founder. For "subrepeat"
cases we take the first `motif_bp` of that founder, repeat it
`copies` times in-place, and keep the rest of the founder unchanged.
The resulting founder is then tiled across ~30 kb and corrupted
with per-base substitution noise.

Outputs:
  <out>.fasta — concatenated FASTA, one record per case
  <out>.truth.tsv — per-record ground truth + injection parameters

Usage:
    python3 tools/subrepeat_sim/inject.py test_data/subrepeat_sim
"""
import random
import sys
from pathlib import Path
from dataclasses import dataclass

# ---------- Case definitions ----------------------------------------------

@dataclass
class Case:
    name: str
    base: str            # "TR" | "HOR3" | "TR_JITTER"
    founder_bp: int      # for TR + TR_JITTER: founder size before injection;
                         # for HOR3: monomer size before injection
    motif_bp: int        # 0 = no subrepeat
    copies: int          # number of subrepeat copies injected (=1 means
                         # no actual tandem — the motif sits once with the
                         # original founder tail beyond)
    expected: str        # "false" | "true" | "borderline" | "degenerate"


CASES = [
    # Negative controls -----------------------------------------------------
    Case("C01_tr_clean",         "TR",         180,  0, 0, "false"),
    Case("C02_hor_clean",        "HOR3",       100,  0, 0, "false"),
    Case("C10_tr_clean_v2",      "TR",         180,  0, 0, "false"),
    Case("C11_tr_jitter",        "TR_JITTER", 2000,  0, 0, "false"),
    # Subrepeat injected with increasing occupancy -------------------------
    Case("C03_tr_sub36_x1",      "TR",         180, 36, 1, "borderline"),
    Case("C04_tr_sub36_x2",      "TR",         180, 36, 2, "true"),
    Case("C05_tr_sub36_x3",      "TR",         180, 36, 3, "true"),
    Case("C06_tr_sub36_x4",      "TR",         180, 36, 4, "true"),
    Case("C07_tr_sub36_x5",      "TR",         180, 36, 5, "degenerate"),
    # Larger founder + subrepeat -------------------------------------------
    Case("C08_tr_300_sub60_x2",  "TR",         300, 60, 2, "true"),
    # HOR with subrepeat inside each monomer slot --------------------------
    Case("C09_hor_sub30_x2",     "HOR3",       100, 30, 2, "true"),
]


# ---------- Sequence helpers ----------------------------------------------

ALPHABET = b"ACGT"


def random_seq(n: int, rng: random.Random) -> bytes:
    return bytes(rng.choice(ALPHABET) for _ in range(n))


def mutate(seq: bytes, sub_rate: float, rng: random.Random) -> bytes:
    out = bytearray(seq)
    for i in range(len(out)):
        if rng.random() < sub_rate:
            out[i] = rng.choice(ALPHABET)
    return bytes(out)


def inject_subrepeat(founder: bytes, motif_bp: int, copies: int,
                     per_copy_sub_rate: float, rng: random.Random) -> bytes:
    """Replace the first `motif_bp · copies` bp of the founder with
    `copies` tandem copies of the founder's leading `motif_bp` bp,
    each copy independently mutated at `per_copy_sub_rate`.
    The remaining `founder_bp − motif_bp · copies` bp (the "gap")
    stays as in the original founder.

    Founder size is preserved — this matches TRC_104's structure
    (founder = subrepeat region + gap, total = original founder).

    Per-copy divergence mirrors real biology: each subrepeat copy
    has accumulated independent substitutions over time. Without
    this, all sampled-pair identities collapse to ~1.0 and the
    bimodal-distribution gates in rescore never fire.

    If `motif_bp · copies ≥ founder_bp`, the founder becomes a pure
    tandem of the motif (occupancy = 1.0; degenerate case).
    """
    if motif_bp == 0 or copies == 0:
        return founder
    master_motif = founder[:motif_bp]
    diverged_copies = b"".join(
        mutate(master_motif, per_copy_sub_rate, rng) for _ in range(copies)
    )
    subrepeat_region = diverged_copies[: len(founder)]
    gap_start = min(len(subrepeat_region), len(founder))
    gap = founder[gap_start:]
    return subrepeat_region + gap


# ---------- Per-base-type builders ----------------------------------------

def build_tr(case: Case, target_len: int, sub_rate: float,
             per_copy_sub_rate: float, jitter_bp: int, seed: int):
    """TR with optional per-copy boundary jitter.

    With `jitter_bp > 0`, every founder copy independently has its
    length perturbed by `random.randint(-jitter_bp, jitter_bp)` — small
    pad / trim at the boundary, modelling real-array indel-driven
    boundary drift. This is essential for measuring metrics that should
    tolerate phase jitter (e.g. the founder-period autocorrelation).
    """
    rng = random.Random(seed)
    base_founder = random_seq(case.founder_bp, rng)
    new_founder = inject_subrepeat(base_founder, case.motif_bp, case.copies,
                                   per_copy_sub_rate, rng)
    nominal_len = len(new_founder)

    if jitter_bp == 0:
        n_copies = max(1, target_len // nominal_len)
        array = new_founder * n_copies
    else:
        pieces = []
        n_copies = 0
        while sum(len(p) for p in pieces) < target_len:
            delta = rng.randint(-jitter_bp, jitter_bp)
            if delta > 0:
                pieces.append(new_founder + random_seq(delta, rng))
            elif delta < 0:
                pieces.append(new_founder[:delta])
            else:
                pieces.append(new_founder)
            n_copies += 1
        array = b"".join(pieces)

    array = mutate(array, sub_rate, rng)
    return array, nominal_len, n_copies, nominal_len  # tile = founder


def build_hor3(case: Case, target_len: int, sub_rate_intra: float,
               sub_rate_inter: float, per_copy_sub_rate: float, seed: int):
    rng = random.Random(seed)
    template = random_seq(case.founder_bp, rng)
    slots = []
    for _ in range(3):
        s = mutate(template, sub_rate_inter, rng)
        s = inject_subrepeat(s, case.motif_bp, case.copies,
                             per_copy_sub_rate, rng)
        slots.append(s)
    slot_len = len(slots[0])
    tile = b"".join(slots)
    tile_len = len(tile)
    n_tiles = max(1, target_len // tile_len)
    array = tile * n_tiles
    array = mutate(array, sub_rate_intra, rng)
    return array, slot_len, n_tiles, tile_len


def build_tr_jitter(case: Case, target_len: int, sub_rate: float,
                    jitter: int, seed: int):
    """TR with per-copy ±jitter bp length variation. Generates kite
    peaks at and near the founder period — useful as a near-founder
    false-positive probe."""
    rng = random.Random(seed)
    base_founder = random_seq(case.founder_bp, rng)
    pieces = []
    n_pieces = 0
    while sum(len(p) for p in pieces) < target_len:
        delta = rng.randint(-jitter, jitter)
        if delta > 0:
            pad = random_seq(delta, rng)
            piece = base_founder + pad
        elif delta < 0:
            piece = base_founder[:delta]
        else:
            piece = base_founder
        pieces.append(piece)
        n_pieces += 1
    array = b"".join(pieces)
    array = mutate(array, sub_rate, rng)
    return array, case.founder_bp, n_pieces, case.founder_bp


# ---------- Main ----------------------------------------------------------

def main():
    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
        "test_data/subrepeat_sim")
    out_dir.mkdir(parents=True, exist_ok=True)

    fa_path = out_dir / "subrepeat_sim.fasta"
    truth_path = out_dir / "truth.tsv"

    TARGET_LEN = 30_000   # ~30 kb
    SUB_RATE = 0.02              # array-wide noise (per-base, post-tiling)
    PER_COPY_SUB_RATE = 0.12     # per-subrepeat-copy divergence inside the founder
                                  # (matches real TRC_104 cov≈0.285 more closely)
    HOR_INTER_RATE = 0.05
    TR_JITTER_BP = 5             # ±5 bp boundary drift per founder copy (realistic)
    BIG_JITTER_BP = 20           # extreme jitter for the C11 near-founder probe
    BASE_SEED = 42

    records = []
    truth_rows = []

    for i, case in enumerate(CASES):
        seed = BASE_SEED + i * 100

        if case.base == "TR":
            seq, founder, n_in_array, tile_bp = build_tr(
                case, TARGET_LEN, SUB_RATE, PER_COPY_SUB_RATE,
                TR_JITTER_BP, seed)
        elif case.base == "HOR3":
            seq, founder, n_in_array, tile_bp = build_hor3(
                case, TARGET_LEN, SUB_RATE, HOR_INTER_RATE,
                PER_COPY_SUB_RATE, seed)
        elif case.base == "TR_JITTER":
            seq, founder, n_in_array, tile_bp = build_tr_jitter(
                case, TARGET_LEN, SUB_RATE, BIG_JITTER_BP, seed)
        else:
            raise ValueError(f"Unknown base: {case.base}")

        records.append((case.name, seq))

        occupancy = ((case.motif_bp * case.copies) / founder
                     if case.motif_bp > 0 else 0.0)
        truth_rows.append({
            "record_id": case.name,
            "base_type": case.base,
            "expected_subrepeat": case.expected,
            "true_founder_bp": founder,
            "true_tile_bp": tile_bp,
            "n_tile_in_array": n_in_array,
            "true_subrepeat_motif_bp": case.motif_bp if case.motif_bp > 0 else "NA",
            "true_subrepeat_copies": case.copies if case.motif_bp > 0 else 0,
            "subrepeat_occupancy_frac": round(occupancy, 4),
            "array_length_bp": len(seq),
        })

    # Write FASTA (80-char wrapped).
    with open(fa_path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for j in range(0, len(seq), 80):
                f.write(seq[j:j + 80].decode("ascii") + "\n")

    # Write truth TSV.
    headers = list(truth_rows[0].keys())
    with open(truth_path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for r in truth_rows:
            f.write("\t".join(str(r[h]) for h in headers) + "\n")

    total_bp = sum(len(s) for _, s in records)
    print(f"Wrote {fa_path} — {len(records)} records, {total_bp} bp total")
    print(f"Wrote {truth_path}")


if __name__ == "__main__":
    main()
