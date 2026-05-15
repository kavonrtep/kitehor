#!/usr/bin/env python3
"""Wider HOR test grid (~1600 cases), test3.

Goes broader than test2 in three directions:

  1. Wider monomer-length range          (15 - 2000 bp, vs 30 - 1500)
  2. Larger HOR multiplicities           (up to N = 12, vs 8)
  3. Two new strata that target known    null_submono, hor_with_submono
     algorithm failure modes              very_short, large_N

Strata
------
  null              350  N=1 plain TR, full monomer/drift range
  null_homog        100  N=1 + heavy monomer conversion
  null_submono      200  N=1 with monomer built from sub-motif x k
                         (the CEN6 TRC_3 pattern; must NOT be called HOR)
  hor_clean         250  N>=2, no conversions, broad
  hor_strong        200  N>=2 + block-conversion fraction 0.2 - 1.5
  hor_decayed       150  N>=2 + monomer conversion
  boundary          100  N>=2 with very weak intra divergence
  stress             50  high indels + extreme monomer / array sizes
  hor_with_submono  100  N>=2 where the base monomer ALSO has internal
                         sub-periodicity (does the detector still call
                         HOR correctly when sub-monomer signal is present?)
  very_short         50  monomer 15-50 bp, includes SSR-adjacent cases
  large_N            50  N in {8, 10, 12}

Total: 1600 cases.

submono_k column
----------------
The new ``submono_k`` parameter triggers the sub-motif-tiling path in
simulate_hor.simulate_array: when submono_k >= 2, the base monomer is
built by tiling a random motif of length ``monomer_len // submono_k``
``submono_k`` times. The autocorrelation spectrum then has a real peak
at ``monomer_len / submono_k`` even though there is no biological HOR.
"""

import csv
import math
import random
from pathlib import Path

OUT = Path(__file__).parent / "params.tsv"
RNG = random.Random(20260512)

FIELDS = [
    "case_id", "monomer_len", "hor_order", "n_blocks",
    "sub_rate_intra", "sub_rate_inter",
    "indel_rate_intra", "indel_rate_inter",
    "block_conversions", "monomer_conversions",
    "submono_k", "seed",
]


def loguniform(lo, hi, rng):
    return math.exp(rng.uniform(math.log(lo), math.log(hi)))


def n_blocks_for(mlen, hor, kb_lo, kb_hi, rng):
    target_bp = loguniform(kb_lo * 1000, kb_hi * 1000, rng)
    return max(2, round(target_bp / (mlen * hor)))


# ---- per-stratum generators ----

def gen_null(rng):
    mlen = round(loguniform(15, 2000, rng))
    return {
        "monomer_len": mlen, "hor_order": 1,
        "n_blocks": n_blocks_for(mlen, 1, 15, 100, rng),
        "sub_rate_intra": 0.0,
        "sub_rate_inter": rng.uniform(0.003, 0.10),
        "indel_rate_intra": 0.0,
        "indel_rate_inter": rng.uniform(0.0, 0.025),
        "block_conv_frac": 0.0, "monomer_conv_frac": 0.0, "submono_k": 1,
    }


def gen_null_homog(rng):
    d = gen_null(rng)
    d["monomer_conv_frac"] = rng.uniform(0.5, 2.0)
    return d


def gen_null_submono(rng):
    """N=1 plain TR but with internal sub-periodicity inside the monomer."""
    k = rng.choice([2, 3, 4, 5])
    sub_len = round(loguniform(15, 200, rng))
    mlen = sub_len * k
    return {
        "monomer_len": mlen, "hor_order": 1,
        "n_blocks": n_blocks_for(mlen, 1, 15, 100, rng),
        "sub_rate_intra": 0.0,
        "sub_rate_inter": rng.uniform(0.005, 0.08),
        "indel_rate_intra": 0.0,
        "indel_rate_inter": rng.uniform(0.0, 0.02),
        "block_conv_frac": 0.0, "monomer_conv_frac": 0.0, "submono_k": k,
    }


def gen_hor_clean(rng):
    mlen = round(loguniform(30, 2000, rng))
    hor = rng.choice([2, 3, 4, 5, 6, 8, 10])
    return {
        "monomer_len": mlen, "hor_order": hor,
        "n_blocks": n_blocks_for(mlen, hor, 15, 100, rng),
        "sub_rate_intra": rng.uniform(0.02, 0.25),
        "sub_rate_inter": rng.uniform(0.003, 0.08),
        "indel_rate_intra": rng.uniform(0.0, 0.02),
        "indel_rate_inter": rng.uniform(0.0, 0.02),
        "block_conv_frac": 0.0, "monomer_conv_frac": 0.0, "submono_k": 1,
    }


def gen_hor_strong(rng):
    d = gen_hor_clean(rng)
    d["block_conv_frac"] = rng.uniform(0.2, 1.5)
    return d


def gen_hor_decayed(rng):
    d = gen_hor_clean(rng)
    d["monomer_conv_frac"] = rng.uniform(0.1, 0.7)
    if rng.random() < 0.5:
        d["block_conv_frac"] = rng.uniform(0.1, 0.5)
    return d


def gen_boundary(rng):
    mlen = round(loguniform(50, 500, rng))
    hor = rng.choice([2, 3, 4, 5, 6])
    return {
        "monomer_len": mlen, "hor_order": hor,
        "n_blocks": n_blocks_for(mlen, hor, 20, 60, rng),
        "sub_rate_intra": rng.uniform(0.005, 0.03),
        "sub_rate_inter": rng.uniform(0.01, 0.05),
        "indel_rate_intra": rng.uniform(0.0, 0.015),
        "indel_rate_inter": rng.uniform(0.0, 0.015),
        "block_conv_frac": 0.0, "monomer_conv_frac": 0.0, "submono_k": 1,
    }


def gen_stress(rng):
    mlen = round(loguniform(30, 1500, rng))
    is_hor = rng.random() < 0.5
    hor = rng.choice([2, 3, 4, 6]) if is_hor else 1
    return {
        "monomer_len": mlen, "hor_order": hor,
        "n_blocks": n_blocks_for(mlen, hor, 8, 50, rng),
        "sub_rate_intra": rng.uniform(0.05, 0.15) if is_hor else 0.0,
        "sub_rate_inter": rng.uniform(0.02, 0.08),
        "indel_rate_intra": rng.uniform(0.02, 0.05) if is_hor else 0.0,
        "indel_rate_inter": rng.uniform(0.02, 0.05),
        "block_conv_frac": 0.0, "monomer_conv_frac": 0.0, "submono_k": 1,
    }


def gen_hor_with_submono(rng):
    """HOR-N where the BASE monomer also has internal sub-periodicity.
    Should still be called as HOR-N (the dominant signal is the HOR
    period, not the sub-monomer)."""
    sub_k = rng.choice([2, 3, 4])
    sub_len = round(loguniform(30, 200, rng))
    mlen = sub_len * sub_k
    hor = rng.choice([2, 3, 4, 5, 6])
    return {
        "monomer_len": mlen, "hor_order": hor,
        "n_blocks": n_blocks_for(mlen, hor, 20, 80, rng),
        "sub_rate_intra": rng.uniform(0.05, 0.20),
        "sub_rate_inter": rng.uniform(0.005, 0.05),
        "indel_rate_intra": rng.uniform(0.0, 0.015),
        "indel_rate_inter": rng.uniform(0.0, 0.015),
        "block_conv_frac": 0.0, "monomer_conv_frac": 0.0,
        "submono_k": sub_k,
    }


def gen_very_short(rng):
    mlen = rng.randint(15, 50)
    is_hor = rng.random() < 0.5
    hor = rng.choice([2, 3, 4]) if is_hor else 1
    return {
        "monomer_len": mlen, "hor_order": hor,
        "n_blocks": n_blocks_for(mlen, hor, 10, 50, rng),
        "sub_rate_intra": rng.uniform(0.02, 0.20) if is_hor else 0.0,
        "sub_rate_inter": rng.uniform(0.005, 0.05),
        "indel_rate_intra": rng.uniform(0.0, 0.01) if is_hor else 0.0,
        "indel_rate_inter": rng.uniform(0.0, 0.01),
        "block_conv_frac": 0.0, "monomer_conv_frac": 0.0, "submono_k": 1,
    }


def gen_large_N(rng):
    mlen = round(loguniform(80, 500, rng))
    hor = rng.choice([8, 10, 12])
    return {
        "monomer_len": mlen, "hor_order": hor,
        "n_blocks": n_blocks_for(mlen, hor, 20, 80, rng),
        "sub_rate_intra": rng.uniform(0.05, 0.20),
        "sub_rate_inter": rng.uniform(0.005, 0.05),
        "indel_rate_intra": rng.uniform(0.0, 0.015),
        "indel_rate_inter": rng.uniform(0.0, 0.015),
        "block_conv_frac": 0.0, "monomer_conv_frac": 0.0, "submono_k": 1,
    }


STRATA = [
    ("null",          350, gen_null),
    ("nullhomog",     100, gen_null_homog),
    ("nullsubmono",   200, gen_null_submono),
    ("horclean",      250, gen_hor_clean),
    ("horstrong",     200, gen_hor_strong),
    ("hordecayed",    150, gen_hor_decayed),
    ("boundary",      100, gen_boundary),
    ("stress",         50, gen_stress),
    ("horsubmono",    100, gen_hor_with_submono),
    ("veryshort",      50, gen_very_short),
    ("largeN",         50, gen_large_N),
]


def main():
    rows = []
    for label, n, gen in STRATA:
        for i in range(n):
            d = gen(RNG)
            n_blocks = d["n_blocks"]
            n_monomers = n_blocks * d["hor_order"]
            rows.append({
                "case_id": f"{label}_{i:04d}",
                "monomer_len": d["monomer_len"],
                "hor_order": d["hor_order"],
                "n_blocks": n_blocks,
                "sub_rate_intra": round(d["sub_rate_intra"], 5),
                "sub_rate_inter": round(d["sub_rate_inter"], 5),
                "indel_rate_intra": round(d["indel_rate_intra"], 5),
                "indel_rate_inter": round(d["indel_rate_inter"], 5),
                "block_conversions": round(n_blocks * d["block_conv_frac"]),
                "monomer_conversions": round(n_monomers * d["monomer_conv_frac"]),
                "submono_k": d["submono_k"],
                "seed": "",
            })

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {len(rows)} cases to {OUT}")


if __name__ == "__main__":
    main()
