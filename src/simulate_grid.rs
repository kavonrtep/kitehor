//! Grid orchestrator: read a params.tsv, simulate every case, write
//! sequences.fasta, truth.tsv, monomers.tsv, events.tsv, alternatives.tsv.
//!
//! Schema matches `ground_truth/simulate_hor.py` for the four legacy
//! files; `alternatives.tsv` is new and surfaces alternative-valid
//! (tile, multiplicity, founder) hierarchies for cases where the
//! synthetic data has more than one biologically valid description
//! (`submono_k >= 2` cases).

use crate::errors::Result;
use crate::simulate::{simulate, SimulateParams};
use anyhow::Context;
use rayon::prelude::*;
use std::io::Write;
use std::path::Path;

/// One row of the input params.tsv.
#[derive(Debug, Clone)]
pub struct GridRow {
    pub case_id: String,
    pub monomer_len: usize,
    pub hor_order: usize,
    pub n_blocks: usize,
    pub sub_rate_intra: f64,
    pub sub_rate_inter: f64,
    pub indel_rate_intra: f64,
    pub indel_rate_inter: f64,
    pub block_conversions: usize,
    pub monomer_conversions: usize,
    pub submono_k: usize,
    pub seed_override: Option<u64>,
}

/// An alternative-valid hierarchy: a (tile, multiplicity, founder)
/// triple that the data also legitimately exhibits, alongside the
/// primary truth.
#[derive(Debug, Clone)]
pub struct Alternative {
    pub rank: usize,
    pub tile: usize,
    pub multiplicity: usize,
    pub founder: usize,
    pub kind: &'static str,
}

/// Enumerate alternative-valid hierarchies for a single case.
///
/// Rules (deterministic, based on `submono_k` and `hor_order`):
///
///   submono_k = 1, hor_order = 1   →  primary only
///                  (tile = monomer_len, k=1, founder = monomer_len)
///
///   submono_k = 1, hor_order >= 2  →  primary only
///                  (tile = monomer_len * hor_order, k = hor_order,
///                   founder = monomer_len)
///
///   submono_k >= 2, hor_order = 1  →  3 rank-ordered alternatives:
///       rank 1: tile = monomer_len, k = 1, founder = monomer_len
///       rank 2: tile = monomer_len, k = submono_k, founder = sub_len
///       rank 3: tile = sub_len,     k = 1,         founder = sub_len
///
///   submono_k >= 2, hor_order >= 2 →  3 rank-ordered alternatives:
///       rank 1: tile = monomer_len * hor_order, k = hor_order,
///                  founder = monomer_len            (primary)
///       rank 2: tile = monomer_len * hor_order,
///                  k = hor_order * submono_k,
///                  founder = sub_len                (sub-motif level)
///       rank 3: tile = monomer_len, k = submono_k,
///                  founder = sub_len                (sub-motif HOR)
///
/// Predictions matching any rank are scored as "valid", though rank-1
/// (the simulator's primary truth) is the preferred answer.
pub fn alternatives_for(row: &GridRow) -> Vec<Alternative> {
    let primary_tile = row.monomer_len * row.hor_order;
    let primary = Alternative {
        rank: 1,
        tile: primary_tile,
        multiplicity: row.hor_order,
        founder: row.monomer_len,
        kind: "primary",
    };
    if row.submono_k < 2 {
        return vec![primary];
    }
    let sub_len = row.monomer_len / row.submono_k;
    if row.hor_order == 1 {
        // nullsubmono-style: a TR with a sub-motif tiled inside.
        vec![
            primary,
            Alternative {
                rank: 2,
                tile: row.monomer_len,
                multiplicity: row.submono_k,
                founder: sub_len,
                kind: "submotif_hor_inside_monomer",
            },
            Alternative {
                rank: 3,
                tile: sub_len,
                multiplicity: 1,
                founder: sub_len,
                kind: "submotif_as_monomer",
            },
        ]
    } else {
        // horsubmono-style: HOR-N where the founder also has internal
        // sub-period.
        vec![
            primary,
            Alternative {
                rank: 2,
                tile: primary_tile,
                multiplicity: row.hor_order * row.submono_k,
                founder: sub_len,
                kind: "submotif_level_inside_horunit",
            },
            Alternative {
                rank: 3,
                tile: row.monomer_len,
                multiplicity: row.submono_k,
                founder: sub_len,
                kind: "submotif_hor_inside_monomer",
            },
        ]
    }
}

/// FNV-1a 64-bit hash, used to derive a deterministic per-case seed
/// when the row's `seed` column is blank.
fn derive_seed(master: u64, case_id: &str) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in master.to_le_bytes().iter() {
        h ^= *b as u64;
        h = h.wrapping_mul(0x100_0000_01b3);
    }
    h ^= b':' as u64;
    h = h.wrapping_mul(0x100_0000_01b3);
    for b in case_id.as_bytes() {
        h ^= *b as u64;
        h = h.wrapping_mul(0x100_0000_01b3);
    }
    h | 1 // Xorshift requires non-zero state.
}

/// Parse params.tsv into a list of GridRow.
pub fn parse_params(path: &Path) -> Result<Vec<GridRow>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .map_err(|e| crate::errors::HordetectError::InvalidParam(format!("{e}")))?;
    let headers: Vec<String> = rdr
        .headers()
        .map_err(|e| crate::errors::HordetectError::InvalidParam(format!("{e}")))?
        .iter()
        .map(|s| s.to_string())
        .collect();
    let idx = |name: &str| -> Result<usize> {
        headers.iter().position(|h| h == name).ok_or_else(|| {
            crate::errors::HordetectError::InvalidParam(format!(
                "params.tsv missing column `{name}`"
            ))
        })
    };
    let i_case = idx("case_id")?;
    let i_ml = idx("monomer_len")?;
    let i_ho = idx("hor_order")?;
    let i_nb = idx("n_blocks")?;
    let i_si = idx("sub_rate_intra")?;
    let i_se = idx("sub_rate_inter")?;
    let i_ii = idx("indel_rate_intra")?;
    let i_ie = idx("indel_rate_inter")?;
    let i_bc = idx("block_conversions")?;
    let i_mc = idx("monomer_conversions")?;
    let i_sk = idx("submono_k")?;
    let i_seed = headers.iter().position(|h| h == "seed");

    let mut rows = Vec::new();
    for (lineno, rec) in rdr.records().enumerate() {
        let rec = rec.map_err(|e| {
            crate::errors::HordetectError::InvalidParam(format!("line {}: {e}", lineno + 2))
        })?;
        let parse_int = |i: usize| -> Result<usize> {
            rec.get(i)
                .ok_or_else(|| {
                    crate::errors::HordetectError::InvalidParam(format!(
                        "line {}: missing column index {i}",
                        lineno + 2
                    ))
                })?
                .trim()
                .parse::<usize>()
                .map_err(|e| {
                    crate::errors::HordetectError::InvalidParam(format!(
                        "line {}: int parse error: {e}",
                        lineno + 2
                    ))
                })
        };
        let parse_f64 = |i: usize| -> Result<f64> {
            rec.get(i)
                .ok_or_else(|| {
                    crate::errors::HordetectError::InvalidParam(format!(
                        "line {}: missing column index {i}",
                        lineno + 2
                    ))
                })?
                .trim()
                .parse::<f64>()
                .map_err(|e| {
                    crate::errors::HordetectError::InvalidParam(format!(
                        "line {}: float parse error: {e}",
                        lineno + 2
                    ))
                })
        };
        let case_id = rec
            .get(i_case)
            .map(|s| s.trim().to_string())
            .unwrap_or_default();
        if case_id.is_empty() {
            continue;
        }
        let seed_override = i_seed.and_then(|i| {
            rec.get(i).and_then(|s| {
                let t = s.trim();
                if t.is_empty() {
                    None
                } else {
                    t.parse::<u64>().ok()
                }
            })
        });
        rows.push(GridRow {
            case_id,
            monomer_len: parse_int(i_ml)?,
            hor_order: parse_int(i_ho)?,
            n_blocks: parse_int(i_nb)?,
            sub_rate_intra: parse_f64(i_si)?,
            sub_rate_inter: parse_f64(i_se)?,
            indel_rate_intra: parse_f64(i_ii)?,
            indel_rate_inter: parse_f64(i_ie)?,
            block_conversions: parse_int(i_bc)?,
            monomer_conversions: parse_int(i_mc)?,
            submono_k: parse_int(i_sk).unwrap_or(1).max(1),
            seed_override,
        });
    }
    Ok(rows)
}

/// Per-case bundle returned from the parallel simulate stage. Strings
/// are pre-formatted so the serial writer just appends.
struct CaseOutput {
    fasta_record: String,
    truth_row: String,
    monomer_rows: String,
    event_rows: String,
    alt_rows: String,
}

const TRUTH_HEADER: &str =
    "case_id\tmonomer_len\thor_order\tn_blocks\tsub_rate_intra\tsub_rate_inter\t\
     indel_rate_intra\tindel_rate_inter\tblock_conversions\tmonomer_conversions\t\
     submono_k\tseed\tarray_length\tn_monomers\t\
     mean_intra_block_id\tmean_homologous_id\tmean_cross_position_id\thor_signal";

const MONOMERS_HEADER: &str = "case_id\tmonomer_idx\tblock_idx\tfounder_idx\tstart\tend\tlength";

const EVENTS_HEADER: &str = "case_id\tevent_order\tscope\tsource_idx\ttarget_idx";

const ALTERNATIVES_HEADER: &str = "case_id\trank\ttile\tmultiplicity\tfounder\tkind";

fn fmt_f64(v: f64) -> String {
    if v.is_nan() {
        "NA".to_string()
    } else {
        format!("{v:.4}")
    }
}

fn simulate_one(row: &GridRow, master_seed: u64) -> Result<CaseOutput> {
    let seed = row
        .seed_override
        .unwrap_or_else(|| derive_seed(master_seed, &row.case_id));
    let params = SimulateParams {
        monomer_len: row.monomer_len,
        hor_order: row.hor_order,
        n_blocks: row.n_blocks,
        sub_rate_intra: row.sub_rate_intra,
        sub_rate_inter: row.sub_rate_inter,
        indel_rate_intra: row.indel_rate_intra,
        indel_rate_inter: row.indel_rate_inter,
        block_conversions: row.block_conversions,
        monomer_conversions: row.monomer_conversions,
        submono_k: row.submono_k,
        seed,
    };
    let (array, truth, monomers, events) = simulate(&row.case_id, &params)?;

    // FASTA record (60-char wrap).
    let mut fasta = String::with_capacity(array.length + array.length / 60 + 64);
    fasta.push('>');
    fasta.push_str(&array.id);
    fasta.push('\n');
    for chunk in array.seq.chunks(60) {
        // SAFETY: array.seq is uppercase ACGT (plus N) — valid UTF-8 bytes.
        fasta.push_str(std::str::from_utf8(chunk).unwrap_or(""));
        fasta.push('\n');
    }

    let truth_row = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        truth.case_id,
        truth.monomer_len,
        truth.hor_order,
        truth.n_blocks,
        fmt_f64(truth.sub_rate_intra),
        fmt_f64(truth.sub_rate_inter),
        fmt_f64(truth.indel_rate_intra),
        fmt_f64(truth.indel_rate_inter),
        truth.block_conversions,
        truth.monomer_conversions,
        truth.submono_k,
        truth.seed,
        truth.array_length,
        truth.n_monomers,
        fmt_f64(truth.mean_intra_block_id),
        fmt_f64(truth.mean_homologous_id),
        fmt_f64(truth.mean_cross_position_id),
        fmt_f64(truth.hor_signal),
    );

    let mut monomer_rows = String::with_capacity(monomers.len() * 40);
    for (idx, m) in monomers.iter().enumerate() {
        monomer_rows.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            row.case_id,
            idx,
            m.block_idx,
            m.founder_idx,
            m.start,
            m.end,
            m.end - m.start,
        ));
    }

    let mut event_rows = String::new();
    for e in &events {
        event_rows.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\n",
            row.case_id, e.event_order, e.scope, e.source_idx, e.target_idx
        ));
    }

    let mut alt_rows = String::new();
    for a in alternatives_for(row) {
        alt_rows.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            row.case_id, a.rank, a.tile, a.multiplicity, a.founder, a.kind
        ));
    }

    Ok(CaseOutput {
        fasta_record: fasta,
        truth_row,
        monomer_rows,
        event_rows,
        alt_rows,
    })
}

/// Run the full grid: parse params.tsv, simulate all cases in parallel,
/// write the five output files in `outdir`.
pub fn run_grid(params_path: &Path, outdir: &Path, master_seed: u64) -> anyhow::Result<()> {
    std::fs::create_dir_all(outdir).with_context(|| format!("creating outdir {outdir:?}"))?;
    let rows = parse_params(params_path)
        .map_err(|e| anyhow::anyhow!("{e}"))
        .with_context(|| format!("reading {params_path:?}"))?;

    let outputs: Vec<Result<CaseOutput>> = rows
        .par_iter()
        .map(|row| simulate_one(row, master_seed))
        .collect();

    let fasta_path = outdir.join("sequences.fasta");
    let truth_path = outdir.join("truth.tsv");
    let monomers_path = outdir.join("monomers.tsv");
    let events_path = outdir.join("events.tsv");
    let alt_path = outdir.join("alternatives.tsv");

    let mut fasta_f = std::fs::File::create(&fasta_path)?;
    let mut truth_f = std::fs::File::create(&truth_path)?;
    let mut mono_f = std::fs::File::create(&monomers_path)?;
    let mut ev_f = std::fs::File::create(&events_path)?;
    let mut alt_f = std::fs::File::create(&alt_path)?;

    writeln!(truth_f, "{TRUTH_HEADER}")?;
    writeln!(mono_f, "{MONOMERS_HEADER}")?;
    writeln!(ev_f, "{EVENTS_HEADER}")?;
    writeln!(alt_f, "{ALTERNATIVES_HEADER}")?;

    let mut n_ok = 0usize;
    let mut n_err = 0usize;
    for (row, out) in rows.iter().zip(outputs) {
        match out {
            Ok(o) => {
                fasta_f.write_all(o.fasta_record.as_bytes())?;
                truth_f.write_all(o.truth_row.as_bytes())?;
                mono_f.write_all(o.monomer_rows.as_bytes())?;
                ev_f.write_all(o.event_rows.as_bytes())?;
                alt_f.write_all(o.alt_rows.as_bytes())?;
                n_ok += 1;
            }
            Err(e) => {
                eprintln!("[err] {}: {}", row.case_id, e);
                n_err += 1;
            }
        }
    }
    eprintln!("simulate-grid: wrote {n_ok} cases ({n_err} errors) to {outdir:?}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alternatives_primary_only_for_submono_k_1() {
        let row = GridRow {
            case_id: "x".into(),
            monomer_len: 200,
            hor_order: 3,
            n_blocks: 10,
            sub_rate_intra: 0.0,
            sub_rate_inter: 0.0,
            indel_rate_intra: 0.0,
            indel_rate_inter: 0.0,
            block_conversions: 0,
            monomer_conversions: 0,
            submono_k: 1,
            seed_override: None,
        };
        let alts = alternatives_for(&row);
        assert_eq!(alts.len(), 1);
        assert_eq!(alts[0].rank, 1);
        assert_eq!(alts[0].tile, 600);
        assert_eq!(alts[0].multiplicity, 3);
        assert_eq!(alts[0].founder, 200);
    }

    #[test]
    fn alternatives_for_nullsubmono() {
        // nullsubmono case: monomer_len=200, hor_order=1, submono_k=4.
        // sub_len = 50.
        let row = GridRow {
            case_id: "x".into(),
            monomer_len: 200,
            hor_order: 1,
            n_blocks: 10,
            sub_rate_intra: 0.0,
            sub_rate_inter: 0.0,
            indel_rate_intra: 0.0,
            indel_rate_inter: 0.0,
            block_conversions: 0,
            monomer_conversions: 0,
            submono_k: 4,
            seed_override: None,
        };
        let alts = alternatives_for(&row);
        assert_eq!(alts.len(), 3);
        assert_eq!(alts[0].tile, 200);
        assert_eq!(alts[0].multiplicity, 1);
        assert_eq!(alts[0].founder, 200);
        assert_eq!(alts[1].tile, 200);
        assert_eq!(alts[1].multiplicity, 4);
        assert_eq!(alts[1].founder, 50);
        assert_eq!(alts[2].tile, 50);
        assert_eq!(alts[2].multiplicity, 1);
        assert_eq!(alts[2].founder, 50);
    }

    #[test]
    fn alternatives_for_horsubmono() {
        // horsubmono_0000-style: monomer_len=212, hor_order=3, submono_k=4.
        // sub_len = 53.
        let row = GridRow {
            case_id: "h".into(),
            monomer_len: 212,
            hor_order: 3,
            n_blocks: 54,
            sub_rate_intra: 0.0,
            sub_rate_inter: 0.0,
            indel_rate_intra: 0.0,
            indel_rate_inter: 0.0,
            block_conversions: 0,
            monomer_conversions: 0,
            submono_k: 4,
            seed_override: None,
        };
        let alts = alternatives_for(&row);
        assert_eq!(alts.len(), 3);
        // Primary: HOR-unit at 212*3 = 636.
        assert_eq!(alts[0].tile, 636);
        assert_eq!(alts[0].multiplicity, 3);
        assert_eq!(alts[0].founder, 212);
        // Sub-motif inside HOR-unit: tile 636, k=12, founder=53.
        assert_eq!(alts[1].tile, 636);
        assert_eq!(alts[1].multiplicity, 12);
        assert_eq!(alts[1].founder, 53);
        // Sub-motif HOR inside founder: tile=212, k=4, founder=53.
        assert_eq!(alts[2].tile, 212);
        assert_eq!(alts[2].multiplicity, 4);
        assert_eq!(alts[2].founder, 53);
    }

    #[test]
    fn derive_seed_is_deterministic_and_different_across_cases() {
        let s1 = derive_seed(42, "horclean_0000");
        let s2 = derive_seed(42, "horclean_0000");
        let s3 = derive_seed(42, "horclean_0001");
        assert_eq!(s1, s2);
        assert_ne!(s1, s3);
        assert_ne!(s1, 0);
    }
}
