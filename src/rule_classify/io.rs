//! TSV I/O for `rule-classify`. Reads kite peaks long-format,
//! emits verdicts.tsv with the prototype's 10-column schema and
//! `%.6g` float formatting (byte-equivalent with
//! `tools/rule_proto/rule_proto.py`).

use super::cluster::Cluster;
use super::decide::{PeakRow, Verdict, VerdictKind};
use anyhow::{anyhow, Context, Result};
use std::io::Write;
use std::path::{Path, PathBuf};

/// Output filename for the verdicts TSV (`<prefix>.verdicts.tsv`).
pub fn verdicts_path(out_prefix: &Path) -> PathBuf {
    let mut p = out_prefix.as_os_str().to_owned();
    p.push(".verdicts.tsv");
    PathBuf::from(p)
}

/// Read a kite peaks long-format TSV and group rows by `case_id`,
/// preserving **first-appearance order** (matches pandas
/// `groupby(sort=False)`).
pub fn read_peaks_grouped(path: &Path) -> Result<Vec<(String, Vec<PeakRow>)>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = rdr.headers()?.clone();
    let mut idx = ColumnIndex::default();
    for (i, h) in headers.iter().enumerate() {
        match h {
            "case_id" => idx.case_id = Some(i),
            "rank" => idx.rank = Some(i),
            "period" => idx.period = Some(i),
            "score2_norm" => idx.score2_norm = Some(i),
            _ => {}
        }
    }
    let case_idx = idx
        .case_id
        .ok_or_else(|| anyhow!("missing 'case_id' column in {:?}", path))?;
    let rank_idx = idx
        .rank
        .ok_or_else(|| anyhow!("missing 'rank' column in {:?}", path))?;
    let period_idx = idx
        .period
        .ok_or_else(|| anyhow!("missing 'period' column in {:?}", path))?;
    let score_idx = idx
        .score2_norm
        .ok_or_else(|| anyhow!("missing 'score2_norm' column in {:?}", path))?;

    let mut order: Vec<String> = Vec::new();
    let mut by_id: std::collections::HashMap<String, Vec<PeakRow>> =
        std::collections::HashMap::new();
    for (rec_idx, rec) in rdr.records().enumerate() {
        let rec = rec.with_context(|| format!("reading record {} of {:?}", rec_idx, path))?;
        let case_id = rec.get(case_idx).unwrap_or("").to_string();
        let rank: u32 = rec
            .get(rank_idx)
            .unwrap_or("0")
            .trim()
            .parse()
            .with_context(|| format!("parsing rank at row {}", rec_idx))?;
        let period: usize = rec
            .get(period_idx)
            .unwrap_or("0")
            .trim()
            .parse()
            .with_context(|| format!("parsing period at row {}", rec_idx))?;
        let score2_norm: f64 = rec
            .get(score_idx)
            .unwrap_or("0")
            .trim()
            .parse()
            .with_context(|| format!("parsing score2_norm at row {}", rec_idx))?;
        by_id
            .entry(case_id.clone())
            .or_insert_with(|| {
                order.push(case_id.clone());
                Vec::new()
            })
            .push(PeakRow {
                rank,
                period,
                score2_norm,
            });
    }
    Ok(order
        .into_iter()
        .map(|id| {
            let v = by_id.remove(&id).unwrap_or_default();
            (id, v)
        })
        .collect())
}

#[derive(Default)]
struct ColumnIndex {
    case_id: Option<usize>,
    rank: Option<usize>,
    period: Option<usize>,
    score2_norm: Option<usize>,
}

/// Write the verdicts TSV. Format matches `tools/rule_proto/rule_proto.py`:
/// 10 columns, header on line 1, `%.6g` float format with empty cells
/// for `None`.
pub fn write_verdicts(path: &Path, verdicts: &[Verdict]) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut f = std::fs::File::create(path).with_context(|| format!("creating {:?}", path))?;
    writeln!(
        f,
        "case_id\tverdict\tfounder\tmultiplicity\ttile\tfounder_score\ttile_score\tconfidence\tn_clusters\treason"
    )?;
    for v in verdicts {
        writeln!(
            f,
            "{cid}\t{vk}\t{f_}\t{m}\t{t}\t{fs}\t{ts}\t{c}\t{nc}\t{r}",
            cid = v.case_id,
            vk = match v.kind {
                VerdictKind::Hor => "hor",
                VerdictKind::SimpleTr => "simple_tr",
                VerdictKind::Unresolved => "unresolved",
            },
            f_ = fmt_opt_g(6, v.founder),
            m = fmt_opt_u32(v.multiplicity),
            t = fmt_opt_g(6, v.tile),
            fs = fmt_opt_g(6, v.founder_score),
            ts = fmt_opt_g(6, v.tile_score),
            c = fmt_opt_g(6, v.confidence),
            nc = v.n_clusters,
            r = v.reason,
        )?;
    }
    Ok(())
}

/// Write a `<case_id>.clusters.tsv` to `dir`. Columns:
/// `rep_period total_score score_frac n_peaks min_rank periods`.
pub fn write_clusters_dump(dir: &Path, case_id: &str, clusters: &[Cluster]) -> Result<()> {
    if clusters.is_empty() {
        return Ok(());
    }
    let max_s = clusters
        .iter()
        .map(|c| c.total_score)
        .fold(f64::NEG_INFINITY, f64::max);
    let path = dir.join(format!("{case_id}.clusters.tsv"));
    let mut f = std::fs::File::create(&path).with_context(|| format!("creating {:?}", &path))?;
    writeln!(
        f,
        "rep_period\ttotal_score\tscore_frac\tn_peaks\tmin_rank\tperiods"
    )?;
    let mut sorted: Vec<&Cluster> = clusters.iter().collect();
    sorted.sort_by(|a, b| {
        b.total_score
            .partial_cmp(&a.total_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    for c in sorted {
        let rep = (c.rep_period * 100.0).round() / 100.0;
        let frac = if max_s > 0.0 { c.total_score / max_s } else { 0.0 };
        let periods_csv: String = c
            .periods
            .iter()
            .map(|p| p.to_string())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}",
            fmt_g(6, rep),
            fmt_g(6, c.total_score),
            fmt_g(6, frac),
            c.n_peaks,
            c.min_rank,
            periods_csv,
        )?;
    }
    Ok(())
}

fn fmt_opt_u32(v: Option<u32>) -> String {
    match v {
        Some(x) => x.to_string(),
        None => String::new(),
    }
}

fn fmt_opt_g(precision: usize, v: Option<f64>) -> String {
    match v {
        Some(x) => fmt_g(precision, x),
        None => String::new(),
    }
}

/// Python `%g`-style float formatting. Used by:
/// - rule-classify (`%.6g`)
/// - summary-merge (`%.4g`)
/// - hor-validate (`%.6g`)
///
/// Semantics:
/// - precision = number of significant digits
/// - choose scientific when `e < -4 || e >= precision`, else fixed
/// - strip trailing zeros from the fractional part; drop dangling `.`
/// - scientific uses lowercase `e`, signed exponent, min 2-digit width
pub fn fmt_g(precision: usize, x: f64) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x.is_sign_negative() { "-inf" } else { "inf" }.to_string();
    }
    if x == 0.0 {
        return "0".to_string();
    }
    let abs = x.abs();
    let e = abs.log10().floor() as i32;
    let p = precision as i32;
    if e >= -4 && e < p {
        let dec = (p - 1 - e).max(0) as usize;
        let s = format!("{:.*}", dec, x);
        strip_trailing_zeros(&s)
    } else {
        let dec = (p - 1).max(0) as usize;
        let raw = format!("{:.*e}", dec, x);
        python_style_scientific(&raw)
    }
}

fn strip_trailing_zeros(s: &str) -> String {
    if !s.contains('.') {
        return s.to_string();
    }
    let trimmed = s.trim_end_matches('0');
    trimmed.trim_end_matches('.').to_string()
}

/// Rust `{:.*e}` yields e.g. `1.23e6`, `1.23e-5`. Python emits
/// `1.23e+06`, `1.23e-05`. Reformat: strip trailing zeros from
/// mantissa, sign the exponent, pad exponent to ≥ 2 digits.
fn python_style_scientific(raw: &str) -> String {
    let Some((mantissa, exp_part)) = raw.split_once('e') else {
        return raw.to_string();
    };
    let mantissa_clean = strip_trailing_zeros(mantissa);
    let (sign, num) = if let Some(stripped) = exp_part.strip_prefix('-') {
        ('-', stripped)
    } else if let Some(stripped) = exp_part.strip_prefix('+') {
        ('+', stripped)
    } else {
        ('+', exp_part)
    };
    let n: i64 = num.parse().unwrap_or(0);
    format!("{mantissa_clean}e{sign}{:02}", n)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check_g(p: usize, x: f64, want: &str) {
        let got = fmt_g(p, x);
        assert_eq!(got, want, "fmt_g({p}, {x}) = {got} (expected {want})");
    }

    #[test]
    fn fmt_g_basics() {
        check_g(6, 0.0, "0");
        check_g(6, 1.0, "1");
        check_g(6, 100.0, "100");
        check_g(6, 300.12, "300.12");
        check_g(6, 0.5, "0.5");
        // 6 sig digits: round
        check_g(6, 0.1234567, "0.123457");
        // Scientific lower bound (abs < 1e-4)
        check_g(6, 0.0000123, "1.23e-05");
        // Scientific upper bound (abs >= 10^precision)
        check_g(6, 1234567.0, "1.23457e+06");
    }

    #[test]
    fn fmt_g_precision_4() {
        check_g(4, 0.21333870, "0.2133");
        check_g(4, 1234.0, "1234");
        check_g(4, 12345.0, "1.234e+04");
    }

    #[test]
    fn fmt_g_negatives() {
        check_g(6, -300.12, "-300.12");
        check_g(6, -0.0000123, "-1.23e-05");
    }

    #[test]
    fn fmt_g_special() {
        assert_eq!(fmt_g(6, f64::NAN), "nan");
        assert_eq!(fmt_g(6, f64::INFINITY), "inf");
        assert_eq!(fmt_g(6, f64::NEG_INFINITY), "-inf");
    }
}
