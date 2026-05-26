//! TSV merge + writer for `kitehor summary-merge`. Reads three input
//! TSVs (verdicts, tandem_validate, ssr) as `HashMap<column, value>`
//! rows, performs an outer join on `record_id`/`case_id`, applies
//! per-column formatting rules to match pandas's
//! `to_csv(float_format="%.4g")`.
//!
//! The merged schema dropped `length_bp` and all
//! subrepeat/hor_validate-derived columns (per the unified-detector
//! port in commit 2 of the tandem_validate plan). Users who need
//! `length_bp` can join on `record_id` against `<prefix>.kite.tsv`.

use super::{combined_class, Config};
use crate::rule_classify::io::fmt_g;
use anyhow::{anyhow, Context, Result};
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Standard target column order. Optional SSR diagnostic columns are
/// filtered out at the end if absent from the source TSV.
pub const TARGET_COLUMNS: &[&str] = &[
    "record_id",
    "hor_verdict",
    "hor_founder",
    "hor_multiplicity",
    "hor_tile",
    "hor_confidence",
    // tandem_validate
    "tv_decision",
    "tv_host_period",
    "tv_best_candidate_period",
    "tv_best_candidate_kind",
    "tv_density",
    "tv_spatial_contrast",
    "tv_phase_contrast",
    "tv_n_windows_total",
    "tv_n_windows_present",
    "tv_reason",
    // ssr (required)
    "ssr_flag",
    "ssr_dominant_motif",
    "ssr_dominant_motif_length",
    "ssr_dominant_motif_repeats",
    "ssr_dominant_motif_coverage_pct",
    "ssr_total_coverage_pct",
    "ssr_top_motifs",
    // ssr diagnostic / consensus columns (optional)
    "ssr_method",
    "consensus_period_bp",
    "consensus_monomer",
    "ssr_raw_dominant_motif",
    "ssr_raw_dominant_motif_coverage_pct",
    "ssr_raw_total_coverage_pct",
    "ssr_raw_n_regions",
    "ssr_raw_top_motifs",
    "combined_class",
];

/// Per-column formatting hint for the writer. Mirrors pandas's dtype
/// inference + `float_format="%.4g"` semantics.
#[derive(Debug, Clone, Copy)]
enum ColFmt {
    /// Pass through (object dtype: strings, mixed numeric/NA). Empty
    /// for missing.
    Passthrough,
    /// Pass through; if absent, write `default`.
    StringDefault(&'static str),
    /// Parse as float, format with `%.4g`. Empty for missing.
    Float4g,
    /// Parse as float, format with `%.4g`. Default `default` when missing.
    Float4gDefault(f64),
    /// Integer column passthrough. Empty for missing.
    Int,
}

fn col_fmt(col: &str) -> ColFmt {
    match col {
        "record_id" => ColFmt::Passthrough,
        "hor_verdict" => ColFmt::StringDefault("unresolved"),
        "hor_founder" => ColFmt::Float4g,
        "hor_multiplicity" => ColFmt::Float4g, // float64 in pandas (NaN for non-HOR/missing)
        "hor_tile" => ColFmt::Float4g,
        "hor_confidence" => ColFmt::Float4g,
        // tandem_validate. `tv_decision` and `tv_reason` are strings;
        // numeric columns are Passthrough so the upstream %.6g rounding
        // from tandem_validate.tsv is preserved verbatim — re-parsing
        // and re-formatting at %.4g would silently drop precision.
        "tv_decision" => ColFmt::StringDefault(""),
        "tv_best_candidate_kind" => ColFmt::Passthrough,
        "tv_reason" => ColFmt::Passthrough,
        "tv_host_period" => ColFmt::Passthrough,
        "tv_best_candidate_period" => ColFmt::Passthrough,
        "tv_density" => ColFmt::Passthrough,
        "tv_spatial_contrast" => ColFmt::Passthrough,
        "tv_phase_contrast" => ColFmt::Passthrough,
        "tv_n_windows_total" => ColFmt::Int,
        "tv_n_windows_present" => ColFmt::Int,
        // ssr (required)
        "ssr_flag" => ColFmt::StringDefault("no"),
        "ssr_dominant_motif" => ColFmt::Passthrough,
        "ssr_dominant_motif_length" => ColFmt::Float4g,
        "ssr_dominant_motif_repeats" => ColFmt::Int,
        "ssr_dominant_motif_coverage_pct" => ColFmt::Float4gDefault(0.0),
        "ssr_total_coverage_pct" => ColFmt::Float4g,
        "ssr_top_motifs" => ColFmt::Passthrough,
        // ssr optional
        "ssr_method" => ColFmt::Passthrough,
        "consensus_period_bp" => ColFmt::Float4g,
        "consensus_monomer" => ColFmt::Passthrough,
        "ssr_raw_dominant_motif" => ColFmt::Passthrough,
        "ssr_raw_dominant_motif_coverage_pct" => ColFmt::Float4g,
        "ssr_raw_total_coverage_pct" => ColFmt::Float4g,
        "ssr_raw_n_regions" => ColFmt::Int,
        "ssr_raw_top_motifs" => ColFmt::Passthrough,
        "combined_class" => ColFmt::Passthrough,
        _ => ColFmt::Passthrough,
    }
}

pub fn summary_path(out_prefix: &Path) -> PathBuf {
    let mut p = out_prefix.as_os_str().to_owned();
    p.push(".summary.tsv");
    PathBuf::from(p)
}

/// Result of merging the inputs.
pub struct Merged {
    /// Columns actually present (subset of [`TARGET_COLUMNS`]).
    pub columns: Vec<&'static str>,
    /// One row per record, keyed by column name. Missing values are
    /// represented by the absence of the key.
    pub rows: Vec<HashMap<String, String>>,
}

/// Outer-join all inputs on record_id. ID order is the union of all
/// three sources sorted lexicographically (matches the prototype's
/// pandas `merge(how="outer")` row order).
pub fn merge_inputs(
    verdicts_path: &Path,
    tandem_validate_path: &Path,
    ssr_path: &Path,
    cfg: &Config,
) -> Result<Merged> {
    let (hor_hdr, hor_rows) = read_tsv(verdicts_path)?;
    let (tv_hdr, tv_rows) = read_tsv(tandem_validate_path)?;
    let (ssr_hdr, ssr_rows) = read_tsv(ssr_path)?;

    let hor_idx = column_index(&hor_hdr, "case_id")
        .ok_or_else(|| anyhow!("verdicts.tsv missing 'case_id'"))?;
    let tv_idx = column_index(&tv_hdr, "record_id")
        .ok_or_else(|| anyhow!("tandem_validate.tsv missing 'record_id'"))?;
    let ssr_idx = column_index(&ssr_hdr, "record_id")
        .ok_or_else(|| anyhow!("ssr.tsv missing 'record_id'"))?;

    let hor_map = build_hor(&hor_hdr, &hor_rows, hor_idx)?;
    let tv_map = build_tv(&tv_hdr, &tv_rows, tv_idx)?;
    let (ssr_map, has_ssr_optional) = build_ssr(&ssr_hdr, &ssr_rows, ssr_idx);

    // ID order: union across the three sources, lexicographically sorted.
    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut id_order: Vec<String> = Vec::new();
    for source in [
        &tv_rows
            .iter()
            .map(|r| r[tv_idx].clone())
            .collect::<Vec<_>>(),
        &hor_rows
            .iter()
            .map(|r| r[hor_idx].clone())
            .collect::<Vec<_>>(),
        &ssr_rows
            .iter()
            .map(|r| r[ssr_idx].clone())
            .collect::<Vec<_>>(),
    ] {
        for id in source {
            if seen.insert(id.clone()) {
                id_order.push(id.clone());
            }
        }
    }
    id_order.sort();

    let mut out_rows: Vec<HashMap<String, String>> = Vec::with_capacity(id_order.len());
    for id in &id_order {
        let mut r: HashMap<String, String> = HashMap::new();
        r.insert("record_id".into(), id.clone());
        for src in [&tv_map, &hor_map, &ssr_map] {
            if let Some(m) = src.get(id) {
                for (k, v) in m {
                    r.insert(k.clone(), v.clone());
                }
            }
        }
        let hor_v = r
            .get("hor_verdict")
            .map(String::as_str)
            .unwrap_or("unresolved");
        // v0.11: cascade reads the array-scale raw total (not the
        // potentially-inflated consensus-path dominant_motif_coverage_pct).
        let raw_total_pct = r
            .get("ssr_raw_total_coverage_pct")
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(0.0);
        let tv_d = r.get("tv_decision").map(String::as_str).unwrap_or("");
        let cls = combined_class(hor_v, raw_total_pct, tv_d, cfg);
        r.insert("combined_class".into(), cls.to_string());
        out_rows.push(r);
    }

    let columns: Vec<&'static str> = TARGET_COLUMNS
        .iter()
        .copied()
        .filter(|c| match *c {
            // Always emitted
            "record_id" | "combined_class" => true,
            "hor_verdict" | "hor_founder" | "hor_multiplicity" | "hor_tile" | "hor_confidence" => {
                true
            }
            "tv_decision"
            | "tv_host_period"
            | "tv_best_candidate_period"
            | "tv_best_candidate_kind"
            | "tv_density"
            | "tv_spatial_contrast"
            | "tv_phase_contrast"
            | "tv_n_windows_total"
            | "tv_n_windows_present"
            | "tv_reason" => true,
            "ssr_flag"
            | "ssr_dominant_motif"
            | "ssr_dominant_motif_length"
            | "ssr_dominant_motif_repeats"
            | "ssr_dominant_motif_coverage_pct"
            | "ssr_total_coverage_pct"
            | "ssr_top_motifs"
            // v0.11: ssr_raw_total_coverage_pct is now load-bearing
            // for the cascade, so it's always emitted.
            | "ssr_raw_total_coverage_pct" => true,
            // Optional SSR diagnostic columns — emit only when the source TSV had them.
            "ssr_method"
            | "consensus_period_bp"
            | "consensus_monomer"
            | "ssr_raw_dominant_motif"
            | "ssr_raw_dominant_motif_coverage_pct"
            | "ssr_raw_n_regions"
            | "ssr_raw_top_motifs" => has_ssr_optional.contains(c),
            _ => false,
        })
        .collect();

    Ok(Merged {
        columns,
        rows: out_rows,
    })
}

/// Apply column transforms for verdicts.tsv. `case_id → record_id`;
/// the 5 columns we keep get renamed with `hor_` prefix (except verdict).
fn build_hor(
    hdr: &[String],
    rows: &[Vec<String>],
    id_idx: usize,
) -> Result<HashMap<String, HashMap<String, String>>> {
    let need = [
        ("verdict", "hor_verdict"),
        ("founder", "hor_founder"),
        ("multiplicity", "hor_multiplicity"),
        ("tile", "hor_tile"),
        ("confidence", "hor_confidence"),
    ];
    let mut col_map: Vec<(usize, &'static str)> = Vec::with_capacity(need.len());
    for (src, dst) in need {
        let i = column_index(hdr, src)
            .ok_or_else(|| anyhow!("verdicts.tsv missing column {:?}", src))?;
        col_map.push((i, dst));
    }
    let mut out: HashMap<String, HashMap<String, String>> = HashMap::new();
    for row in rows {
        let id = row[id_idx].clone();
        let mut m: HashMap<String, String> = HashMap::new();
        for &(i, dst) in &col_map {
            let v = row.get(i).cloned().unwrap_or_default();
            if !v.is_empty() {
                m.insert(dst.to_string(), v);
            }
        }
        out.insert(id, m);
    }
    Ok(out)
}

/// Apply column transforms for tandem_validate.tsv. All kept columns
/// get a `tv_` prefix; `decision_hint` → `tv_decision`.
fn build_tv(
    hdr: &[String],
    rows: &[Vec<String>],
    id_idx: usize,
) -> Result<HashMap<String, HashMap<String, String>>> {
    let need: &[(&str, &'static str)] = &[
        ("host_period", "tv_host_period"),
        ("best_candidate_period", "tv_best_candidate_period"),
        ("best_candidate_kind", "tv_best_candidate_kind"),
        ("density", "tv_density"),
        ("spatial_contrast", "tv_spatial_contrast"),
        ("phase_contrast", "tv_phase_contrast"),
        ("n_windows_total", "tv_n_windows_total"),
        ("n_windows_present", "tv_n_windows_present"),
        ("decision_hint", "tv_decision"),
        ("reason", "tv_reason"),
    ];
    let mut col_map: Vec<(usize, &'static str)> = Vec::new();
    for (src, dst) in need {
        let i = column_index(hdr, src)
            .ok_or_else(|| anyhow!("tandem_validate.tsv missing column {:?}", src))?;
        col_map.push((i, dst));
    }
    let mut out: HashMap<String, HashMap<String, String>> = HashMap::new();
    for row in rows {
        let id = row[id_idx].clone();
        let mut m: HashMap<String, String> = HashMap::new();
        for &(i, dst) in &col_map {
            let v = row.get(i).cloned().unwrap_or_default();
            if !v.is_empty() {
                m.insert(dst.to_string(), v);
            }
        }
        out.insert(id, m);
    }
    Ok(out)
}

fn build_ssr(
    hdr: &[String],
    rows: &[Vec<String>],
    id_idx: usize,
) -> (
    HashMap<String, HashMap<String, String>>,
    std::collections::HashSet<&'static str>,
) {
    let required: &[(&str, &'static str)] = &[
        ("ssr_flag", "ssr_flag"),
        ("dominant_motif", "ssr_dominant_motif"),
        ("dominant_motif_length", "ssr_dominant_motif_length"),
        ("dominant_motif_repeats", "ssr_dominant_motif_repeats"),
        (
            "dominant_motif_coverage_pct",
            "ssr_dominant_motif_coverage_pct",
        ),
        ("total_ssr_coverage_pct", "ssr_total_coverage_pct"),
        ("top_motifs", "ssr_top_motifs"),
    ];
    let optional: &[&'static str] = &[
        "ssr_method",
        "consensus_period_bp",
        "consensus_monomer",
        "ssr_raw_dominant_motif",
        "ssr_raw_dominant_motif_coverage_pct",
        "ssr_raw_total_coverage_pct",
        "ssr_raw_n_regions",
        "ssr_raw_top_motifs",
    ];
    let mut col_map: Vec<(usize, &'static str)> = Vec::new();
    for (src, dst) in required {
        if let Some(i) = column_index(hdr, src) {
            col_map.push((i, dst));
        }
    }
    let mut optional_present: std::collections::HashSet<&'static str> =
        std::collections::HashSet::new();
    for &name in optional {
        if let Some(i) = column_index(hdr, name) {
            col_map.push((i, name));
            optional_present.insert(name);
        }
    }
    let mut out: HashMap<String, HashMap<String, String>> = HashMap::new();
    for row in rows {
        let id = row[id_idx].clone();
        let mut m: HashMap<String, String> = HashMap::new();
        for &(i, dst) in &col_map {
            let v = row.get(i).cloned().unwrap_or_default();
            if !v.is_empty() {
                m.insert(dst.to_string(), v);
            }
        }
        out.insert(id, m);
    }
    (out, optional_present)
}

/// Read a TSV. Returns (header columns, row values). Cell values that
/// pandas's default `na_values` set would parse as missing are converted
/// to the empty string here so downstream formatting renders them as
/// empty cells (matches `to_csv(na_rep="")`).
fn read_tsv(path: &Path) -> Result<(Vec<String>, Vec<Vec<String>>)> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {:?}", path))?;
    let mut lines = text.lines();
    let Some(header_line) = lines.next() else {
        return Ok((Vec::new(), Vec::new()));
    };
    let header: Vec<String> = header_line.split('\t').map(String::from).collect();
    let mut rows: Vec<Vec<String>> = Vec::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let cells: Vec<String> = line
            .split('\t')
            .map(|s| {
                if is_pandas_na(s) {
                    "".to_string()
                } else {
                    s.to_string()
                }
            })
            .collect();
        rows.push(cells);
    }
    Ok((header, rows))
}

fn is_pandas_na(s: &str) -> bool {
    matches!(
        s,
        "" | "#N/A"
            | "#N/A N/A"
            | "#NA"
            | "-1.#IND"
            | "-1.#QNAN"
            | "-NaN"
            | "-nan"
            | "1.#IND"
            | "1.#QNAN"
            | "<NA>"
            | "N/A"
            | "NA"
            | "NULL"
            | "NaN"
            | "None"
            | "n/a"
            | "nan"
            | "null"
    )
}

fn column_index(hdr: &[String], name: &str) -> Option<usize> {
    hdr.iter().position(|h| h == name)
}

/// Write the merged TSV with per-column formatting matching pandas's
/// `to_csv(float_format="%.4g", na_rep="")`.
pub fn write_summary(path: &Path, merged: &Merged) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut f = std::fs::File::create(path).with_context(|| format!("creating {:?}", path))?;
    writeln!(f, "{}", merged.columns.join("\t"))?;
    for row in &merged.rows {
        let line: Vec<String> = merged
            .columns
            .iter()
            .map(|&c| format_cell(c, row.get(c).map(String::as_str)))
            .collect();
        writeln!(f, "{}", line.join("\t"))?;
    }
    Ok(())
}

fn format_cell(col: &str, val: Option<&str>) -> String {
    match col_fmt(col) {
        ColFmt::Passthrough => val.unwrap_or("").to_string(),
        ColFmt::StringDefault(d) => match val {
            Some(s) if !s.is_empty() => s.to_string(),
            _ => d.to_string(),
        },
        ColFmt::Int => match val {
            Some(s) if !s.is_empty() => s.to_string(),
            _ => String::new(),
        },
        ColFmt::Float4g => match val.and_then(|s| s.parse::<f64>().ok()) {
            Some(x) => fmt_g(4, x),
            None => String::new(),
        },
        ColFmt::Float4gDefault(d) => match val.and_then(|s| s.parse::<f64>().ok()) {
            Some(x) => fmt_g(4, x),
            None => fmt_g(4, d),
        },
    }
}
