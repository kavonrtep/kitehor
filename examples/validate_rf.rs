//! Cross-validate the Rust RF + Platt against an external reference TSV
//! (typically `predict_verdict.R` output). Useful when regenerating the
//! shipped model artifacts and you want to confirm the Rust loader is
//! still bit-equivalent to ranger.
//!
//! Usage:
//!   cargo run --release --example validate_rf -- \
//!       --features <path>            (input: extract_features.py output)
//!       --reference <path>           (input: predict_verdict.R output)
//!       [--config <classifier.toml>] (default: baked-in)
//!       [--model <hor_score.json>]   (default: models/hor_score.rftrees.json)

use std::collections::HashMap;
use std::path::PathBuf;

use kitehor::classifier::{ClassifierConfig, RandomForest};

fn parse_arg(args: &[String], key: &str, default: Option<&str>) -> Option<String> {
    let mut it = args.iter().peekable();
    while let Some(a) = it.next() {
        if a == key {
            return it.next().cloned();
        }
    }
    default.map(|s| s.to_string())
}

fn parse_f64(s: &str) -> Option<f64> {
    if s == "NA" || s.is_empty() {
        None
    } else {
        s.parse::<f64>().ok()
    }
}

fn read_tsv(path: &str) -> (Vec<String>, Vec<Vec<String>>) {
    let txt = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read tsv {path:?}: {e}"));
    let mut lines = txt.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .map(|s| s.to_string())
        .collect();
    let rows: Vec<Vec<String>> = lines
        .map(|l| l.split('\t').map(|s| s.to_string()).collect())
        .collect();
    (header, rows)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let features_path =
        parse_arg(&args, "--features", None).expect("--features <path> is required");
    let reference_path =
        parse_arg(&args, "--reference", None).expect("--reference <path> is required");
    let model_path = parse_arg(&args, "--model", Some("models/hor_score.rftrees.json")).unwrap();
    let config_path = parse_arg(&args, "--config", None);

    let cfg = match config_path.as_deref() {
        Some(p) => ClassifierConfig::load(PathBuf::from(p)).expect("load config"),
        None => ClassifierConfig::default_baked().expect("baked config"),
    };
    let h_model = RandomForest::load_json(&model_path).expect("load RF");
    let platt = cfg.platt();

    let (hdr_feat, rows_feat) = read_tsv(&features_path);
    let (hdr_pred, rows_pred) = read_tsv(&reference_path);
    let col_feat: HashMap<&str, usize> = hdr_feat
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();
    let col_pred: HashMap<&str, usize> = hdr_pred
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();

    let case_col_p = col_pred["case_id"];
    let r_raw_col = col_pred["hor_score_raw"];
    let r_cal_col = col_pred["hor_score"];
    let mut pred_by_case: HashMap<&str, &Vec<String>> = HashMap::new();
    for r in &rows_pred {
        pred_by_case.insert(r[case_col_p].as_str(), r);
    }
    let case_col_f = col_feat["case_id"];

    let mut max_raw = 0.0_f64;
    let mut max_cal = 0.0_f64;
    let mut sum_raw = 0.0_f64;
    let mut sum_cal = 0.0_f64;
    let mut n = 0usize;
    let mut shown = 0usize;
    for row in &rows_feat {
        let case = row[case_col_f].as_str();
        let mut x = Vec::with_capacity(h_model.feature_names.len());
        for name in &h_model.feature_names {
            let c = col_feat
                .get(name.as_str())
                .unwrap_or_else(|| panic!("missing feature column: {name}"));
            let mut v = parse_f64(&row[*c]);
            if v.is_none() {
                v = match name.as_str() {
                    "h_d1" => Some(cfg.imputation.h_d1),
                    "h_founder" => Some(cfg.imputation.h_founder),
                    _ => None,
                };
            }
            x.push(v.unwrap_or(f64::NAN));
        }
        let raw = h_model.predict(&x);
        let cal = platt.calibrate(raw);
        let Some(p) = pred_by_case.get(case) else {
            eprintln!("[warn] case missing in reference: {case}");
            continue;
        };
        let r_raw: f64 = p[r_raw_col].parse().unwrap_or(f64::NAN);
        let r_cal: f64 = p[r_cal_col].parse().unwrap_or(f64::NAN);
        let d_raw = (raw - r_raw).abs();
        let d_cal = (cal - r_cal).abs();
        max_raw = max_raw.max(d_raw);
        max_cal = max_cal.max(d_cal);
        sum_raw += d_raw;
        sum_cal += d_cal;
        n += 1;
        if shown < 5 || d_raw > 1e-3 {
            println!(
                "{case:<55} R_raw={r_raw:.10} Rust_raw={raw:.10} Δ={d_raw:.2e}   R_cal={r_cal:.10} Rust_cal={cal:.10} Δ={d_cal:.2e}"
            );
            shown += 1;
        }
    }
    println!(
        "\n=== {n} records ===\nmax |Δ raw| = {max_raw:.3e}\nmean |Δ raw|= {:.3e}\nmax |Δ cal| = {max_cal:.3e}\nmean |Δ cal|= {:.3e}",
        sum_raw / n as f64,
        sum_cal / n as f64,
    );
}
