//! TSV merge + writer for `kitehor summary-merge`. Reads four input
//! TSVs as `HashMap<column, value>` rows, performs an outer join on
//! `record_id`/`case_id`, applies per-column formatting rules to match
//! pandas's `to_csv(float_format="%.4g")` byte-for-byte.

use super::{combined_class, Config};
use crate::rule_classify::io::fmt_g;
use anyhow::{anyhow, Context, Result};
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Standard target column order. Optional columns are filtered out at
/// the end if absent from every input source.
pub const TARGET_COLUMNS: &[&str] = &[
    "record_id",
    "length_bp",
    "hor_verdict",
    "hor_founder",
    "hor_multiplicity",
    "hor_tile",
    "hor_confidence",
    "subrepeat_flag",
    "subrepeat_host_period_bp",
    "subrepeat_period_bp",
    "n_subrepeat_blocks",
    "subrepeat_coverage_pct",
    "ssr_flag",
    "ssr_dominant_motif",
    "ssr_dominant_motif_length",
    "ssr_dominant_motif_repeats",
    "ssr_dominant_motif_coverage_pct",
    "ssr_total_coverage_pct",
    "ssr_top_motifs",
    // Diagnostic / consensus columns (optional in ssr.tsv)
    "ssr_method",
    "consensus_period_bp",
    "consensus_monomer",
    "ssr_raw_dominant_motif",
    "ssr_raw_dominant_motif_coverage_pct",
    "ssr_raw_total_coverage_pct",
    "ssr_raw_n_regions",
    "ssr_raw_top_motifs",
    // within-tile validation (optional)
    "density_hint",
    "founder_density",
    "phase_contrast",
    "density_n_windows",
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
        // length_bp: int from subrepeat.tsv. After outer merge missing
        // records gain NaN → pandas coerces to float64 → %.4g. But
        // since length_bp is always populated in normal pipelines,
        // pass through as-is (matches both code paths when no NaN).
        "length_bp" => ColFmt::Int,
        "hor_verdict" => ColFmt::StringDefault("unresolved"),
        "hor_founder" => ColFmt::Float4g,
        "hor_multiplicity" => ColFmt::Float4g, // float64 in pandas (NaN for non-HOR/missing)
        "hor_tile" => ColFmt::Float4g,
        "hor_confidence" => ColFmt::Float4g,
        "subrepeat_flag" => ColFmt::StringDefault("none"),
        // pandas reads "NA" as NaN → column becomes float64 → %.4g applies.
        "subrepeat_host_period_bp" => ColFmt::Float4g,
        "subrepeat_period_bp" => ColFmt::Float4g,
        "n_subrepeat_blocks" => ColFmt::Int,
        "subrepeat_coverage_pct" => ColFmt::Float4g,
        "ssr_flag" => ColFmt::StringDefault("no"),
        "ssr_dominant_motif" => ColFmt::Passthrough,
        "ssr_dominant_motif_length" => ColFmt::Float4g,
        "ssr_dominant_motif_repeats" => ColFmt::Int,
        "ssr_dominant_motif_coverage_pct" => ColFmt::Float4gDefault(0.0),
        "ssr_total_coverage_pct" => ColFmt::Float4g,
        "ssr_top_motifs" => ColFmt::Passthrough,
        "ssr_method" => ColFmt::Passthrough,
        "consensus_period_bp" => ColFmt::Float4g,
        "consensus_monomer" => ColFmt::Passthrough,
        "ssr_raw_dominant_motif" => ColFmt::Passthrough,
        "ssr_raw_dominant_motif_coverage_pct" => ColFmt::Float4g,
        "ssr_raw_total_coverage_pct" => ColFmt::Float4g,
        "ssr_raw_n_regions" => ColFmt::Int,
        "ssr_raw_top_motifs" => ColFmt::Passthrough,
        "density_hint" => ColFmt::StringDefault(""),
        "founder_density" => ColFmt::Passthrough, // "NA" or %.6g float
        "phase_contrast" => ColFmt::Passthrough,
        "density_n_windows" => ColFmt::Int,
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
    /// Columns actually present (subset of `TARGET_COLUMNS`).
    pub columns: Vec<&'static str>,
    /// One row per record, keyed by column name. Missing values are
    /// represented by the absence of the key.
    pub rows: Vec<HashMap<String, String>>,
}

/// Outer-join all inputs on record_id. Walk order:
/// 1. Records in subrepeat.tsv (in file order).
/// 2. Records only in hor (verdicts.tsv) — added at end.
/// 3. Records only in ssr.tsv — added at end after that.
/// 4. within-tile records get joined into existing rows (left join);
///    no new rows added.
pub fn merge_inputs(
    verdicts_path: &Path,
    subrepeat_path: &Path,
    ssr_path: &Path,
    within_tile_path: Option<&Path>,
    cfg: &Config,
) -> Result<Merged> {
    // 1. Read all 4 TSVs as (header, rows).
    let (hor_hdr, hor_rows) = read_tsv(verdicts_path)?;
    let (sub_hdr, sub_rows) = read_tsv(subrepeat_path)?;
    let (ssr_hdr, ssr_rows) = read_tsv(ssr_path)?;
    let within_tuple = match within_tile_path {
        Some(p) => Some(read_tsv(p)?),
        None => None,
    };

    // 2. Per-source column transforms:
    //    - hor: case_id → record_id; subset/rename verdict→hor_verdict, etc.
    //    - sub: keep subset, rename host_period_bp → subrepeat_host_period_bp.
    //    - ssr: rename per the prototype's rename map.
    //    - within: keep density_hint + 3 optional cols; left-join only.
    let hor_idx = column_index(&hor_hdr, "case_id")
        .ok_or_else(|| anyhow!("verdicts.tsv missing 'case_id'"))?;
    let sub_idx = column_index(&sub_hdr, "record_id")
        .ok_or_else(|| anyhow!("subrepeat.tsv missing 'record_id'"))?;
    let ssr_idx = column_index(&ssr_hdr, "record_id")
        .ok_or_else(|| anyhow!("ssr.tsv missing 'record_id'"))?;

    // Transformed per-row hash maps keyed by record_id.
    let hor_map = build_hor(&hor_hdr, &hor_rows, hor_idx)?;
    let sub_order = sub_rows
        .iter()
        .map(|r| r[sub_idx].clone())
        .collect::<Vec<String>>();
    let sub_map = build_sub(&sub_hdr, &sub_rows, sub_idx);
    let (ssr_map, has_ssr_optional) = build_ssr(&ssr_hdr, &ssr_rows, ssr_idx);
    let (within_map, has_within_extra) = match within_tuple {
        Some((hdr, rows)) => {
            let idx = column_index(&hdr, "record_id")
                .ok_or_else(|| anyhow!("hor_within_tile.tsv missing 'record_id'"))?;
            build_within(&hdr, &rows, idx)
        }
        None => (HashMap::new(), false),
    };
    let within_supplied = within_tile_path.is_some();
    // density_hint column is always emitted, even when --within-tile is
    // omitted (matches the prototype's `m["density_hint"] = ""` branch).
    let _ = within_supplied;

    // 3. Outer-join walk. Pandas's `pd.merge(..., how="outer")` produces
    //    output in lexicographically-sorted key order in the user's
    //    environment (confirmed empirically against the prototype). Take
    //    the union of all keys and sort.
    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut id_order: Vec<String> = Vec::new();
    for id in &sub_order {
        if seen.insert(id.clone()) {
            id_order.push(id.clone());
        }
    }
    let hor_order = hor_rows.iter().map(|r| r[hor_idx].clone());
    for id in hor_order {
        if seen.insert(id.clone()) {
            id_order.push(id);
        }
    }
    let ssr_order = ssr_rows.iter().map(|r| r[ssr_idx].clone());
    for id in ssr_order {
        if seen.insert(id.clone()) {
            id_order.push(id);
        }
    }
    id_order.sort();

    // 4. Build per-record merged rows.
    let mut out_rows: Vec<HashMap<String, String>> = Vec::with_capacity(id_order.len());
    for id in &id_order {
        let mut r: HashMap<String, String> = HashMap::new();
        r.insert("record_id".into(), id.clone());
        if let Some(m) = sub_map.get(id) {
            for (k, v) in m {
                r.insert(k.clone(), v.clone());
            }
        }
        if let Some(m) = hor_map.get(id) {
            for (k, v) in m {
                r.insert(k.clone(), v.clone());
            }
        }
        if let Some(m) = ssr_map.get(id) {
            for (k, v) in m {
                r.insert(k.clone(), v.clone());
            }
        }
        if let Some(m) = within_map.get(id) {
            for (k, v) in m {
                r.insert(k.clone(), v.clone());
            }
        }
        // Apply defaults for combined_class computation.
        let hor_v = r
            .get("hor_verdict")
            .map(String::as_str)
            .unwrap_or("unresolved");
        let sub_f = r.get("subrepeat_flag").map(String::as_str).unwrap_or("none");
        let ssr_f = r.get("ssr_flag").map(String::as_str).unwrap_or("no");
        let dom_pct = r
            .get("ssr_dominant_motif_coverage_pct")
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(0.0);
        let dens = if within_supplied {
            r.get("density_hint").map(String::as_str).unwrap_or("")
        } else {
            ""
        };
        let cls = combined_class(hor_v, ssr_f, dom_pct, sub_f, dens, cfg);
        r.insert("combined_class".into(), cls.to_string());
        out_rows.push(r);
    }

    // 5. Decide which columns are present. Optional columns are dropped
    //    if no input ever populated them.
    let columns: Vec<&'static str> = TARGET_COLUMNS
        .iter()
        .copied()
        .filter(|c| {
            match *c {
                // Always present from the merge logic
                "record_id" | "length_bp" | "hor_verdict" | "hor_founder"
                | "hor_multiplicity" | "hor_tile" | "hor_confidence" | "subrepeat_flag"
                | "subrepeat_host_period_bp" | "subrepeat_period_bp"
                | "n_subrepeat_blocks" | "subrepeat_coverage_pct" | "ssr_flag"
                | "ssr_dominant_motif" | "ssr_dominant_motif_length"
                | "ssr_dominant_motif_repeats" | "ssr_dominant_motif_coverage_pct"
                | "ssr_total_coverage_pct" | "ssr_top_motifs" => true,
                "combined_class" => true,
                "density_hint" => true,
                "founder_density" | "phase_contrast" | "density_n_windows" => has_within_extra,
                "ssr_method" | "consensus_period_bp" | "consensus_monomer"
                | "ssr_raw_dominant_motif" | "ssr_raw_dominant_motif_coverage_pct"
                | "ssr_raw_total_coverage_pct" | "ssr_raw_n_regions" | "ssr_raw_top_motifs" => {
                    has_ssr_optional.contains(c)
                }
                _ => false,
            }
        })
        .collect();
    Ok(Merged {
        columns,
        rows: out_rows,
    })
}

/// Apply column transforms for verdicts.tsv. case_id → record_id; the
/// 5 columns we keep get renamed with `hor_` prefix (except verdict).
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

fn build_sub(
    hdr: &[String],
    rows: &[Vec<String>],
    id_idx: usize,
) -> HashMap<String, HashMap<String, String>> {
    let mappings: &[(&str, &str)] = &[
        ("length_bp", "length_bp"),
        ("subrepeat_flag", "subrepeat_flag"),
        ("host_period_bp", "subrepeat_host_period_bp"),
        ("subrepeat_period_bp", "subrepeat_period_bp"),
        ("n_subrepeat_blocks", "n_subrepeat_blocks"),
        ("subrepeat_coverage_pct", "subrepeat_coverage_pct"),
    ];
    let mut col_map: Vec<(usize, &'static str)> = Vec::new();
    for (src, dst) in mappings {
        if let Some(i) = column_index(hdr, src) {
            col_map.push((i, dst));
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
    out
}

fn build_ssr(
    hdr: &[String],
    rows: &[Vec<String>],
    id_idx: usize,
) -> (
    HashMap<String, HashMap<String, String>>,
    std::collections::HashSet<&'static str>,
) {
    // Required (renamed)
    let required: &[(&str, &'static str)] = &[
        ("ssr_flag", "ssr_flag"),
        ("dominant_motif", "ssr_dominant_motif"),
        ("dominant_motif_length", "ssr_dominant_motif_length"),
        ("dominant_motif_repeats", "ssr_dominant_motif_repeats"),
        ("dominant_motif_coverage_pct", "ssr_dominant_motif_coverage_pct"),
        ("total_ssr_coverage_pct", "ssr_total_coverage_pct"),
        ("top_motifs", "ssr_top_motifs"),
    ];
    // Optional (no rename)
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

fn build_within(
    hdr: &[String],
    rows: &[Vec<String>],
    id_idx: usize,
) -> (HashMap<String, HashMap<String, String>>, bool) {
    let cols: &[&'static str] = &[
        "density_hint",
        "founder_density",
        "phase_contrast",
        "density_n_windows",
    ];
    let mut col_map: Vec<(usize, &'static str)> = Vec::new();
    let mut has_extra = false;
    for &c in cols {
        if let Some(i) = column_index(hdr, c) {
            col_map.push((i, c));
            if c != "density_hint" {
                has_extra = true;
            }
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
    (out, has_extra)
}

/// Read a TSV. Returns (header columns, row values). Cell values that
/// would be interpreted by `pandas.read_csv` as missing (the default
/// `na_values`: "NA", "NaN", "null", "N/A", "<NA>", "nan", "None", etc.)
/// are converted to the empty string here so downstream formatting
/// renders them as empty cells (matches `to_csv(na_rep="")`).
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
            .map(|s| if is_pandas_na(s) { "".to_string() } else { s.to_string() })
            .collect();
        rows.push(cells);
    }
    Ok((header, rows))
}

/// pandas's default `na_values` set. Without per-column `na_filter=False`
/// these strings are coerced to NaN on read.
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
