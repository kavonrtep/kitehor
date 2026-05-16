//! IO for the line-width detector.
//!
//! - `load_arrays` reads a FASTA via the existing `crate::io::load_fasta`
//!   helper (already normalises to A/C/G/T/N).
//! - `load_periods` reads a `periods.tsv` with the validation rules
//!   from `detect_impl_plan.md §7` (A12) and groups rows by `array_id`.
//! - `write_*` emit the three frozen-at-M0 TSVs.

use crate::detect::types::{
    PeriodCandidate, Properties, Segment, WidthFeatures, PROPERTIES_HEADER,
    SEGMENTS_HEADER, WIDTH_FEATURES_HEADER,
};
use crate::io::{load_fasta, LoadQc, LoadStatus, LoadedRecord};
use crate::sequence::ArrayRecord;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Read a FASTA into normalised `ArrayRecord`s. Records that fail QC
/// are logged but returned anyway so the property writer can mark
/// them `ambiguous`.
pub fn load_arrays(path: &Path) -> Result<Vec<ArrayRecord>> {
    // Use a permissive QC: the detector itself decides per array
    // whether the input is usable. The legacy QC thresholds were
    // tuned for the rule-based classifier on real centromeric data.
    let qc = LoadQc {
        min_array_bp: 0,
        max_n_fraction: 1.0,
    };
    let loaded: Vec<LoadedRecord> = load_fasta(path, qc)
        .with_context(|| format!("loading FASTA {:?}", path))?;
    let mut out = Vec::with_capacity(loaded.len());
    for r in loaded {
        match &r.status {
            LoadStatus::Ok => out.push(r.record),
            _ => {
                log::warn!(
                    "FASTA record {} flagged: {:?} — passing through",
                    r.record.id, r.status
                );
                out.push(r.record);
            }
        }
    }
    Ok(out)
}

/// Read `periods.tsv` and group by `array_id`.
///
/// Multi-record FASTA semantics (A13): if every row's `array_id` is
/// empty, callers pair by file stem (single-record convention). The
/// detector errors out at join time if a FASTA record has no matching
/// period rows.
pub fn load_periods(path: &Path) -> Result<HashMap<String, Vec<PeriodCandidate>>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening periods TSV {:?}", path))?;
    let headers = rdr.headers()?.clone();
    let idx = |col: &str| {
        headers.iter().position(|h| h == col).ok_or_else(|| {
            anyhow::anyhow!("periods TSV {:?} missing required column `{}`", path, col)
        })
    };
    let i_period = idx("period_bp")?;
    let i_score = idx("period_score")?;
    let i_array = headers.iter().position(|h| h == "array_id");
    let i_source = headers.iter().position(|h| h == "source");

    let mut grouped: HashMap<String, Vec<PeriodCandidate>> = HashMap::new();
    let mut dup_check: HashMap<(String, usize), f64> = HashMap::new();

    for (row_idx, rec) in rdr.records().enumerate() {
        let rec = rec.with_context(|| format!("reading row {} of {:?}", row_idx, path))?;
        let array_id = match i_array {
            Some(i) => rec.get(i).unwrap_or("").trim().to_string(),
            None => String::new(),
        };
        let period_bp: usize = rec
            .get(i_period)
            .unwrap_or("")
            .parse()
            .with_context(|| {
                format!("row {} period_bp not an integer in {:?}", row_idx, path)
            })?;
        let period_score: f64 = rec.get(i_score).unwrap_or("0").parse().unwrap_or(0.0);
        let source = i_source
            .map(|i| rec.get(i).unwrap_or("").to_string())
            .unwrap_or_default();

        // Duplicate detection: keep the max score per (array_id, period_bp).
        let key = (array_id.clone(), period_bp);
        match dup_check.get(&key).copied() {
            Some(prev) if prev >= period_score => {
                log::warn!(
                    "duplicate period row dropped: array={:?} period={} (score {} <= {})",
                    array_id, period_bp, period_score, prev
                );
                continue;
            }
            Some(_) => {
                // Higher-score duplicate replaces the earlier entry. Walk the
                // already-stored vec and drop the lower-score copy.
                if let Some(v) = grouped.get_mut(&array_id) {
                    v.retain(|c| c.period_bp != period_bp);
                }
            }
            None => {}
        }
        dup_check.insert(key, period_score);
        grouped.entry(array_id.clone()).or_default().push(PeriodCandidate {
            array_id,
            period_bp,
            period_score,
            source,
        });
    }
    Ok(grouped)
}

/// Join FASTA records with their period candidates. Returns
/// `(array, periods)` pairs. The `default_array_id` is used when the
/// periods TSV had no `array_id` column (single-record convention).
///
/// When `allow_missing` is `false` (default), a FASTA record with no
/// matching period rows is a hard error — see review finding DH5.
/// Pass `--allow-missing-periods` on the CLI to downgrade it to a
/// warning (the array then runs with no candidate widths and ends up
/// `ambiguous`).
pub fn join_arrays_with_periods(
    arrays: Vec<ArrayRecord>,
    mut periods: HashMap<String, Vec<PeriodCandidate>>,
    default_array_id: Option<&str>,
    allow_missing: bool,
) -> Result<Vec<(ArrayRecord, Vec<PeriodCandidate>)>> {
    let mut out = Vec::with_capacity(arrays.len());
    for arr in arrays {
        // Look up by record id first; fall back to default-stem lookup if
        // the periods TSV was single-record (no array_id column).
        let pers = if let Some(v) = periods.remove(&arr.id) {
            v
        } else if let Some(stem) = default_array_id {
            periods.remove(stem).unwrap_or_default()
        } else {
            periods.remove("").unwrap_or_default()
        };
        if pers.is_empty() {
            if allow_missing {
                log::warn!(
                    "no period candidates for FASTA record `{}` — detector will fall back to UnsupportedWidth (`--allow-missing-periods` is set)",
                    arr.id
                );
            } else {
                anyhow::bail!(
                    "no period candidates for FASTA record `{}`; pass `--allow-missing-periods` to downgrade to a warning",
                    arr.id
                );
            }
        }
        out.push((arr, pers));
    }
    Ok(out)
}

// ---------- Output writers ----------

pub fn properties_path(prefix: &Path) -> PathBuf {
    with_ext(prefix, "properties.tsv")
}
pub fn segments_path(prefix: &Path) -> PathBuf {
    with_ext(prefix, "segments.tsv")
}
pub fn width_features_path(prefix: &Path) -> PathBuf {
    with_ext(prefix, "width_features.tsv")
}

pub fn diagnostics_path(prefix: &Path) -> PathBuf {
    with_ext(prefix, "diagnostics.json")
}

fn with_ext(prefix: &Path, ext: &str) -> PathBuf {
    let mut s = prefix.as_os_str().to_owned();
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

pub fn write_properties(prefix: &Path, rows: &[Properties]) -> Result<()> {
    let path = properties_path(prefix);
    ensure_parent(&path)?;
    let mut f = std::fs::File::create(&path)
        .with_context(|| format!("creating {:?}", path))?;
    writeln!(f, "{}", PROPERTIES_HEADER)?;
    for r in rows {
        writeln!(f, "{}", properties_to_tsv(r))?;
    }
    Ok(())
}

pub fn write_segments(prefix: &Path, rows: &[Segment]) -> Result<()> {
    let path = segments_path(prefix);
    ensure_parent(&path)?;
    let mut f = std::fs::File::create(&path)
        .with_context(|| format!("creating {:?}", path))?;
    writeln!(f, "{}", SEGMENTS_HEADER)?;
    for r in rows {
        writeln!(f, "{}", segment_to_tsv(r))?;
    }
    Ok(())
}

pub fn write_width_features(prefix: &Path, rows: &[WidthFeatures]) -> Result<()> {
    let path = width_features_path(prefix);
    ensure_parent(&path)?;
    let mut f = std::fs::File::create(&path)
        .with_context(|| format!("creating {:?}", path))?;
    writeln!(f, "{}", WIDTH_FEATURES_HEADER)?;
    for r in rows {
        writeln!(f, "{}", width_features_to_tsv(r))?;
    }
    Ok(())
}

/// Write a structured per-run summary (`PREFIX.diagnostics.json`)
/// containing the final decision per array, the chosen base width,
/// the supporting evidence, and the top tested widths. DH12.
pub fn write_diagnostics(
    prefix: &Path,
    properties: &[Properties],
    width_features: &[WidthFeatures],
    segments: &[Segment],
) -> Result<()> {
    use serde_json::json;
    let path = diagnostics_path(prefix);
    ensure_parent(&path)?;

    // Group width_features by array_id for compact per-array
    // diagnostic blocks.
    let mut by_array: std::collections::HashMap<&str, Vec<&WidthFeatures>> =
        std::collections::HashMap::new();
    for w in width_features {
        by_array.entry(w.array_id.as_str()).or_default().push(w);
    }
    let mut segs_by_array: std::collections::HashMap<&str, Vec<&Segment>> =
        std::collections::HashMap::new();
    for s in segments {
        segs_by_array.entry(s.array_id.as_str()).or_default().push(s);
    }

    let arrays: Vec<serde_json::Value> = properties
        .iter()
        .map(|p| {
            let widths: Vec<serde_json::Value> = by_array
                .get(p.array_id.as_str())
                .into_iter()
                .flat_map(|v| v.iter())
                .map(|w| {
                    json!({
                        "width_bp": w.width_bp,
                        "rows": w.rows,
                        "column_ic": w.column_ic,
                        "fraction_conserved_columns": w.fraction_conserved_columns,
                        "row_lag1_similarity": w.row_lag1_similarity,
                        "best_lag": w.best_lag,
                        "best_lag_score": w.best_lag_score,
                        "phase_separation": w.phase_separation,
                        "vertical_edge_rate": w.vertical_edge_rate,
                        "mean_shift_bp": w.mean_shift_bp,
                        "wobble_amplitude_bp": w.wobble_amplitude_bp,
                        "n_phase_shifts": w.n_phase_shifts,
                        "irregularity_score": w.irregularity_score,
                        "class_hint": w.class_hint.as_str(),
                    })
                })
                .collect();
            let segs: Vec<serde_json::Value> = segs_by_array
                .get(p.array_id.as_str())
                .into_iter()
                .flat_map(|v| v.iter())
                .map(|s| {
                    json!({
                        "segment_id": s.segment_id,
                        "start_bp": s.start_bp,
                        "end_bp": s.end_bp,
                        "class": s.class.as_str(),
                        "base_width_bp": s.base_width_bp,
                        "hor_k": s.hor_k,
                    })
                })
                .collect();
            json!({
                "array_id": p.array_id,
                "length_bp": p.length_bp,
                "class": p.class.as_str(),
                "base_width_bp": p.base_width_bp,
                "hor_k": p.hor_k,
                "hor_length_bp": p.hor_length_bp,
                "n_complete_copies": p.n_complete_copies,
                "column_conservation": p.column_conservation,
                "phase_separation": p.phase_separation,
                "wobble_amplitude_bp": p.wobble_amplitude_bp,
                "n_phase_shifts": p.n_phase_shifts,
                "phase_shift_positions": p.phase_shift_positions,
                "phase_shift_offsets": p.phase_shift_offsets,
                "irregularity_score": p.irregularity_score,
                "inter_monomer_identity": p.inter_monomer_identity,
                "confidence": p.confidence,
                "reason": p.reason,
                "width_features": widths,
                "segments": segs,
            })
        })
        .collect();

    let doc = json!({
        "schema_version": 1,
        "n_arrays": properties.len(),
        "arrays": arrays,
    });
    let s = serde_json::to_string_pretty(&doc)?;
    std::fs::write(&path, s).with_context(|| format!("creating {:?}", path))?;
    Ok(())
}

fn ensure_parent(p: &Path) -> Result<()> {
    if let Some(parent) = p.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    Ok(())
}

const NA: &str = "NA";

fn opt_usize(o: &Option<usize>) -> String {
    o.map(|x| x.to_string()).unwrap_or_else(|| NA.to_string())
}
fn opt_f64(o: &Option<f64>) -> String {
    o.map(|x| format!("{:.4}", x))
        .unwrap_or_else(|| NA.to_string())
}
fn list_usize(v: &[usize]) -> String {
    if v.is_empty() {
        NA.to_string()
    } else {
        v.iter().map(usize::to_string).collect::<Vec<_>>().join(",")
    }
}
fn list_i64(v: &[i64]) -> String {
    if v.is_empty() {
        NA.to_string()
    } else {
        v.iter().map(i64::to_string).collect::<Vec<_>>().join(",")
    }
}

fn properties_to_tsv(r: &Properties) -> String {
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        r.array_id,
        r.length_bp,
        r.class.as_str(),
        opt_usize(&r.base_width_bp),
        opt_usize(&r.hor_k),
        opt_usize(&r.hor_length_bp),
        opt_usize(&r.n_complete_copies),
        opt_f64(&r.column_conservation),
        opt_f64(&r.phase_separation),
        opt_f64(&r.mean_shift_bp),
        opt_f64(&r.wobble_amplitude_bp),
        opt_f64(&r.wobble_periodicity_bp),
        r.n_phase_shifts,
        list_usize(&r.phase_shift_positions),
        list_i64(&r.phase_shift_offsets),
        opt_f64(&r.irregularity_score),
        opt_f64(&r.inter_monomer_identity),
        opt_f64(&r.confidence),
        r.n_segments,
        // Tab-strip the reason so it doesn't break the TSV.
        r.reason.replace('\t', " "),
    )
}

fn segment_to_tsv(s: &Segment) -> String {
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        s.array_id,
        s.segment_id,
        s.start_bp,
        s.end_bp,
        s.class.as_str(),
        opt_usize(&s.base_width_bp),
        opt_usize(&s.hor_k),
        opt_f64(&s.column_conservation),
        opt_f64(&s.phase_separation),
        opt_f64(&s.wobble_amplitude_bp),
        opt_f64(&s.irregularity_score),
    )
}

fn width_features_to_tsv(w: &WidthFeatures) -> String {
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        w.array_id,
        w.width_bp,
        w.rows,
        opt_f64(&w.column_ic),
        opt_f64(&w.fraction_conserved_columns),
        opt_f64(&w.row_lag1_similarity),
        opt_usize(&w.best_lag),
        opt_f64(&w.best_lag_score),
        opt_f64(&w.phase_separation),
        opt_f64(&w.vertical_edge_rate),
        opt_usize(&w.column_edge_autocorr_k),
        opt_f64(&w.column_edge_autocorr_score),
        opt_f64(&w.mean_shift_bp),
        opt_f64(&w.wobble_amplitude_bp),
        w.n_phase_shifts,
        opt_f64(&w.irregularity_score),
        w.class_hint.as_str(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::detect::types::{Class, ClassHint};
    use std::io::Write;

    #[test]
    fn empty_properties_round_trip_header_then_zero_rows() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("t");
        write_properties(&prefix, &[]).unwrap();
        let s = std::fs::read_to_string(properties_path(&prefix)).unwrap();
        assert_eq!(s.trim_end(), PROPERTIES_HEADER);
    }

    #[test]
    fn placeholder_writes_all_na() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("t");
        let row = Properties::placeholder("arr1", 12345);
        write_properties(&prefix, &[row]).unwrap();
        let s = std::fs::read_to_string(properties_path(&prefix)).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines.len(), 2);
        let fields: Vec<&str> = lines[1].split('\t').collect();
        let header_fields: Vec<&str> = PROPERTIES_HEADER.split('\t').collect();
        assert_eq!(fields.len(), header_fields.len(), "column count must match");
        assert_eq!(fields[0], "arr1");
        assert_eq!(fields[1], "12345");
        assert_eq!(fields[2], "ambiguous");
        // base_width_bp..confidence are all NA
        for i in 3..18 {
            assert!(
                fields[i] == "NA" || fields[i] == "0",
                "field {i} expected NA/0, got `{}`",
                fields[i]
            );
        }
    }

    #[test]
    fn periods_loader_groups_by_array_id() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("periods.tsv");
        std::fs::File::create(&p).unwrap().write_all(
            b"array_id\tperiod_bp\tperiod_score\tsource\n\
              a1\t171\t0.94\ttrue_base\n\
              a1\t2052\t0.88\ttrue_hor_unit\n\
              a2\t170\t0.94\ttrue_base\n",
        ).unwrap();
        let m = load_periods(&p).unwrap();
        assert_eq!(m.get("a1").unwrap().len(), 2);
        assert_eq!(m.get("a2").unwrap().len(), 1);
    }

    #[test]
    fn periods_loader_dedups_keeping_max_score() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("periods.tsv");
        std::fs::File::create(&p).unwrap().write_all(
            b"array_id\tperiod_bp\tperiod_score\tsource\n\
              a1\t171\t0.50\ttrue_base\n\
              a1\t171\t0.94\ttrue_base\n",
        ).unwrap();
        let m = load_periods(&p).unwrap();
        let v = m.get("a1").unwrap();
        assert_eq!(v.len(), 1);
        assert!((v[0].period_score - 0.94).abs() < 1e-9);
    }

    #[test]
    fn periods_loader_rejects_missing_required_column() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("bad.tsv");
        std::fs::File::create(&p).unwrap().write_all(
            b"array_id\tperiod_bp\n\
              a1\t171\n",
        ).unwrap();
        assert!(load_periods(&p).is_err());
    }

    fn dummy_widths() -> Vec<WidthFeatures> {
        vec![WidthFeatures {
            array_id: "arr".into(),
            width_bp: 171,
            rows: 100,
            column_ic: None,
            fraction_conserved_columns: None,
            row_lag1_similarity: None,
            best_lag: None,
            best_lag_score: None,
            phase_separation: None,
            vertical_edge_rate: None,
            column_edge_autocorr_k: None,
            column_edge_autocorr_score: None,
            mean_shift_bp: None,
            wobble_amplitude_bp: None,
            n_phase_shifts: 0,
            irregularity_score: None,
            class_hint: ClassHint::UnsupportedWidth,
        }]
    }

    #[test]
    fn width_features_writer_emits_one_row_per_width() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("t");
        write_width_features(&prefix, &dummy_widths()).unwrap();
        let s = std::fs::read_to_string(width_features_path(&prefix)).unwrap();
        assert!(s.starts_with(WIDTH_FEATURES_HEADER));
        assert_eq!(s.lines().count(), 2);
    }

    #[test]
    fn segments_writer_emits_header_only_when_empty() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("t");
        write_segments(&prefix, &[]).unwrap();
        let s = std::fs::read_to_string(segments_path(&prefix)).unwrap();
        assert_eq!(s.trim_end(), SEGMENTS_HEADER);
    }

    fn dummy_segments() -> Vec<Segment> {
        vec![Segment {
            array_id: "arr".into(),
            segment_id: 1,
            start_bp: 0,
            end_bp: 10000,
            class: Class::Ambiguous,
            base_width_bp: None,
            hor_k: None,
            column_conservation: None,
            phase_separation: None,
            wobble_amplitude_bp: None,
            irregularity_score: None,
        }]
    }

    #[test]
    fn segments_writer_one_row() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("t");
        write_segments(&prefix, &dummy_segments()).unwrap();
        let s = std::fs::read_to_string(segments_path(&prefix)).unwrap();
        assert_eq!(s.lines().count(), 2);
    }
}
