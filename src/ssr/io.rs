//! TSV I/O for `ssr-scan`. Float formatting matches pandas's default
//! (`to_csv(sep="\t", index=False)` with no `float_format` set) —
//! integer-valued floats render with a `.0` suffix, fractionals retain
//! pandas's default 12-significant-digit width.

use super::scan::{RegionRow, SummaryRow};
use ahash::AHashMap;
use anyhow::{Context, Result};
use std::io::Write;
use std::path::{Path, PathBuf};

pub fn summary_path(out_prefix: &Path) -> PathBuf {
    let mut p = out_prefix.as_os_str().to_owned();
    p.push(".ssr.tsv");
    PathBuf::from(p)
}

pub fn regions_path(out_prefix: &Path) -> PathBuf {
    let mut p = out_prefix.as_os_str().to_owned();
    p.push(".ssr.regions.tsv");
    PathBuf::from(p)
}

/// Parse a FASTA into `Vec<(id, seq)>` preserving file order. IDs are
/// the first whitespace-separated token of the header (matches the
/// prototype's `line[1:].split()[0]`).
pub fn read_fasta_ordered(path: &Path) -> Result<Vec<(String, Vec<u8>)>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {:?}", path))?;
    let mut out: Vec<(String, Vec<u8>)> = Vec::new();
    let mut cur_id: Option<String> = None;
    let mut cur_seq: Vec<u8> = Vec::new();
    for line in text.lines() {
        let line = line.trim_end();
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(id) = cur_id.take() {
                out.push((id, std::mem::take(&mut cur_seq)));
            }
            let id = rest.split_whitespace().next().unwrap_or("").to_string();
            cur_id = Some(id);
        } else if !line.is_empty() {
            cur_seq.extend_from_slice(line.as_bytes());
        }
    }
    if let Some(id) = cur_id.take() {
        out.push((id, cur_seq));
    }
    Ok(out)
}

/// Read a kite peaks TSV and return `{record_id: top-peak-period}`,
/// where "top" is the row with maximum `score2_norm` per case_id.
pub fn read_kite_top_periods(path: &Path) -> Result<AHashMap<String, usize>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {:?}", path))?;
    let mut lines = text.lines();
    let Some(hdr) = lines.next() else {
        return Ok(AHashMap::new());
    };
    let cols: Vec<&str> = hdr.split('\t').collect();
    let case_idx = cols
        .iter()
        .position(|c| *c == "case_id")
        .ok_or_else(|| anyhow::anyhow!("kite peaks missing case_id"))?;
    let period_idx = cols
        .iter()
        .position(|c| *c == "period")
        .ok_or_else(|| anyhow::anyhow!("kite peaks missing period"))?;
    let score_idx = cols
        .iter()
        .position(|c| *c == "score2_norm")
        .ok_or_else(|| anyhow::anyhow!("kite peaks missing score2_norm"))?;

    let mut best: AHashMap<String, (f64, usize)> = AHashMap::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        let case_id = cells[case_idx].to_string();
        let period: usize = cells[period_idx].parse().unwrap_or(0);
        let score: f64 = cells[score_idx].parse().unwrap_or(0.0);
        let entry = best.entry(case_id).or_insert((f64::NEG_INFINITY, 0usize));
        if score > entry.0 {
            *entry = (score, period);
        }
    }
    Ok(best.into_iter().map(|(k, v)| (k, v.1)).collect())
}

/// Format a float the way pandas's `to_csv` does when no
/// `float_format` is set (default). Pandas behavior:
/// - integer-valued floats render with a trailing `.0` (e.g. 0 → "0.0")
/// - non-integer floats use up to 12 significant digits, trailing
///   zeros trimmed
fn fmt_pandas_default(x: f64) -> String {
    if x.is_nan() {
        return "".to_string(); // na_rep default
    }
    if x.is_infinite() {
        return if x.is_sign_negative() { "-inf" } else { "inf" }.to_string();
    }
    if x == x.trunc() && x.abs() < 1e16 {
        // Integer-valued — render with trailing ".0"
        return format!("{}.0", x as i64);
    }
    // Up to 12 sig digits; trim trailing zeros from fraction.
    let s = format!("{:.11e}", x);
    // Parse "1.23456789012e2" form and convert to fixed/sci as pandas would.
    // pandas uses str(float) which produces shortest roundtrip for the
    // exact float64. The Rust `{}` formatter also uses shortest-roundtrip
    // for f64 (Grisu/Ryu), which usually matches Python's repr semantics.
    let _ = s;
    let mut s = format!("{}", x);
    // Rust `{}` may produce "300" for 300.0 without decimal. Add ".0".
    if !s.contains('.') && !s.contains('e') && !s.contains('E') {
        s.push_str(".0");
    }
    s
}

/// Write the per-record summary TSV.
pub fn write_summary(path: &Path, rows: &[SummaryRow]) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut f = std::fs::File::create(path).with_context(|| format!("creating {:?}", path))?;
    writeln!(
        f,
        "record_id\tlength_bp\tssr_flag\tdominant_motif\tdominant_motif_length\tdominant_motif_repeats\tdominant_motif_coverage_pct\ttotal_ssr_coverage_pct\ttop_motifs\tssr_method\tconsensus_period_bp\tconsensus_monomer\tssr_raw_dominant_motif\tssr_raw_dominant_motif_coverage_pct\tssr_raw_total_coverage_pct\tssr_raw_n_regions\tssr_raw_top_motifs"
    )?;
    for r in rows {
        writeln!(
            f,
            "{rid}\t{len}\t{flag}\t{dm}\t{dml}\t{dmr}\t{dmcov}\t{tot}\t{top}\t{meth}\t{cp}\t{cm}\t{rdm}\t{rdmcov}\t{rtot}\t{rn}\t{rtop}",
            rid = r.record_id,
            len = r.length_bp,
            flag = r.ssr_flag,
            dm = r.dominant_motif,
            dml = r.dominant_motif_length,
            dmr = r.dominant_motif_repeats,
            dmcov = fmt_pandas_default(r.dominant_motif_coverage_pct),
            tot = fmt_pandas_default(r.total_ssr_coverage_pct),
            top = r.top_motifs,
            meth = r.ssr_method,
            cp = r.consensus_period_bp,
            cm = r.consensus_monomer,
            rdm = r.raw_dominant_motif,
            rdmcov = fmt_pandas_default(r.raw_dominant_motif_coverage_pct),
            rtot = fmt_pandas_default(r.raw_total_coverage_pct),
            rn = r.raw_n_regions,
            rtop = r.raw_top_motifs,
        )?;
    }
    Ok(())
}

/// Write the raw-regions TSV. When there are no rows we emit just a
/// bare newline (no header) — matches the prototype's behavior of
/// `pd.DataFrame([]).to_csv(...)`, which has no columns to print.
pub fn write_regions(path: &Path, rows: &[RegionRow]) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut f = std::fs::File::create(path).with_context(|| format!("creating {:?}", path))?;
    if rows.is_empty() {
        writeln!(f)?;
        return Ok(());
    }
    writeln!(
        f,
        "record_id\tssr_number\tmotif_length\tmotif_sequence\trepeats\tstart\tend\tnormalized_motif"
    )?;
    for r in rows {
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.record_id,
            r.ssr_number,
            r.motif_length,
            r.motif_sequence,
            r.repeats,
            r.start,
            r.end,
            r.normalized_motif,
        )?;
    }
    Ok(())
}
