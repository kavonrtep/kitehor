//! Smoke test: synthesize three minimal input TSVs (verdicts +
//! tandem_validate + ssr), run `kitehor summary-merge` against them,
//! and verify the v0.11 9-rule combined_class dispatch — every verdict
//! base class (`hor`, `tr`, `unresolved`, `tr_with_subrepeat`) has a
//! `_with_ssr` partner, all driven by the array-scale
//! `raw_total_coverage_pct` (not the consensus-overridden dominant
//! coverage).

use std::path::PathBuf;
use std::process::Command;

fn kitehor_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_kitehor"))
}

fn write_file(path: &std::path::Path, body: &str) {
    std::fs::write(path, body).expect("writing input");
}

#[test]
fn summary_merge_nine_rule_dispatch() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();

    // Six records exercising all nine cascade branches.
    let verdicts = dir.join("v.tsv");
    write_file(
        &verdicts,
        "case_id\tverdict\tfounder\tmultiplicity\ttile\tfounder_score\ttile_score\tconfidence\tn_clusters\treason\n\
         r_hor\thor\t100\t3\t300\t0.3\t0.5\t0.8\t4\ttop_is_multiple_of_founder\n\
         r_hor_ssr\thor\t100\t3\t300\t0.3\t0.5\t0.8\t4\ttop_is_multiple_of_founder\n\
         r_tr\tsimple_tr\t150\t1\t150\t0.5\t0.5\t0.7\t3\tlone_significant_cluster\n\
         r_tr_ssr\tsimple_tr\t150\t1\t150\t0.5\t0.5\t0.7\t3\tlone_significant_cluster\n\
         r_unres\tunresolved\t\t\t\t\t\t\t0\tno_clusters\n\
         r_unres_ssr\tunresolved\t\t\t\t\t\t\t0\tno_clusters\n\
         r_pure_ssr\tsimple_tr\t60\t1\t60\t0.4\t0.4\t0.6\t1\tlone_significant_cluster\n\
         r_sub\thor\t200\t3\t600\t0.2\t0.4\t0.7\t3\ttop_is_multiple_of_founder\n\
         r_sub_ssr\thor\t200\t3\t600\t0.2\t0.4\t0.7\t3\ttop_is_multiple_of_founder\n",
    );

    let tandem = dir.join("tv.tsv");
    // tandem_validate columns: record_id, verdict, host_period, multiplicity,
    // window_bp, n_candidates, candidates, best_candidate_period,
    // best_candidate_kind, density, spatial_contrast, phase_contrast,
    // n_windows_total, n_windows_present, decision_hint, reason
    write_file(
        &tandem,
        "record_id\tverdict\thost_period\tmultiplicity\twindow_bp\tn_candidates\tcandidates\tbest_candidate_period\tbest_candidate_kind\tdensity\tspatial_contrast\tphase_contrast\tn_windows_total\tn_windows_present\tdecision_hint\treason\n\
         r_hor\thor\t300\t3\t300\t1\tfounder/100:d=1.000:sc=0.000:pc=nan:uniform\t100\tfounder\t1\t0\tNA\t10\t10\tconfirms_host\tfounder:uniform\n\
         r_hor_ssr\thor\t300\t3\t300\t1\tfounder/100:d=1.000:sc=0.000:pc=nan:uniform\t100\tfounder\t1\t0\tNA\t10\t10\tconfirms_host\tfounder:uniform\n\
         r_tr\tsimple_tr\t150\t1\t\t0\t\t\t\t\t\tNA\t0\t0\tno_candidates\tno_candidates\n\
         r_tr_ssr\tsimple_tr\t150\t1\t\t0\t\t\t\t\t\tNA\t0\t0\tno_candidates\tno_candidates\n\
         r_unres\tunresolved\t\t1\t\t0\t\t\t\t\t\tNA\t0\t0\tno_host\tno_host\n\
         r_unres_ssr\tunresolved\t\t1\t\t0\t\t\t\t\t\tNA\t0\t0\tno_host\tno_host\n\
         r_pure_ssr\tsimple_tr\t60\t1\t\t0\t\t\t\t\t\tNA\t0\t0\tno_candidates\tno_candidates\n\
         r_sub\thor\t600\t3\t200\t1\tfounder/200:d=0.300:sc=0.500:pc=0.000:localized\t200\tfounder\t0.3\t0.5\t0\t15\t5\tlocalized_subrepeat\tfounder:localized\n\
         r_sub_ssr\thor\t600\t3\t200\t1\tfounder/200:d=0.300:sc=0.500:pc=0.000:localized\t200\tfounder\t0.3\t0.5\t0\t15\t5\tlocalized_subrepeat\tfounder:localized\n",
    );

    // SSR TSV columns: record_id, length_bp, ssr_flag, dominant_motif,
    // dominant_motif_length, dominant_motif_repeats,
    // dominant_motif_coverage_pct, total_ssr_coverage_pct, top_motifs,
    // + optional ssr_raw_total_coverage_pct (v0.11 cascade reads this).
    // Note: per-record `ssr_flag` here is just an input passthrough;
    // the cascade decision uses ssr_raw_total_coverage_pct exclusively.
    let ssr = dir.join("s.tsv");
    write_file(
        &ssr,
        "record_id\tlength_bp\tssr_flag\tdominant_motif\tdominant_motif_length\tdominant_motif_repeats\tdominant_motif_coverage_pct\ttotal_ssr_coverage_pct\ttop_motifs\tssr_raw_total_coverage_pct\n\
         r_hor\t5000\tno\tNA\tNA\t0\t0.0\t0.0\tNA\t0.0\n\
         r_hor_ssr\t5000\tyes\tAT\t2\t100\t40.0\t40.0\tAT:40.0%\t40.0\n\
         r_tr\t6000\tno\tNA\tNA\t0\t0.0\t0.0\tNA\t0.0\n\
         r_tr_ssr\t6000\tyes\tAT\t2\t120\t40.0\t40.0\tAT:40.0%\t40.0\n\
         r_unres\t7000\tno\tNA\tNA\t0\t0.0\t0.0\tNA\t0.0\n\
         r_unres_ssr\t7000\tyes\tAT\t2\t140\t40.0\t40.0\tAT:40.0%\t40.0\n\
         r_pure_ssr\t8000\tyes\tAT\t2\t4000\t90.0\t90.0\tAT:90.0%\t90.0\n\
         r_sub\t9000\tno\tNA\tNA\t0\t0.0\t0.0\tNA\t0.0\n\
         r_sub_ssr\t9000\tyes\tAT\t2\t180\t40.0\t40.0\tAT:40.0%\t40.0\n",
    );

    let out_prefix = dir.join("smry");
    let status = Command::new(kitehor_bin())
        .args(["summary-merge", "--verdicts"])
        .arg(&verdicts)
        .arg("--tandem-validate")
        .arg(&tandem)
        .arg("--ssr")
        .arg(&ssr)
        .arg("-o")
        .arg(&out_prefix)
        .status()
        .unwrap();
    assert!(status.success(), "summary-merge failed");

    let out_path = {
        let mut p = out_prefix.into_os_string();
        p.push(".summary.tsv");
        PathBuf::from(p)
    };
    let body = std::fs::read_to_string(&out_path).unwrap();
    let by_id: std::collections::HashMap<&str, &str> = body
        .lines()
        .skip(1)
        .filter(|l| !l.is_empty())
        .map(|l| {
            let cells: Vec<&str> = l.split('\t').collect();
            // last column is combined_class
            (cells[0], *cells.last().unwrap())
        })
        .collect();
    assert_eq!(by_id.get("r_hor"), Some(&"hor"));
    assert_eq!(by_id.get("r_hor_ssr"), Some(&"hor_with_ssr"));
    assert_eq!(by_id.get("r_tr"), Some(&"tr"));
    assert_eq!(by_id.get("r_tr_ssr"), Some(&"tr_with_ssr"));
    assert_eq!(by_id.get("r_unres"), Some(&"unresolved"));
    assert_eq!(by_id.get("r_unres_ssr"), Some(&"unresolved_with_ssr"));
    assert_eq!(by_id.get("r_pure_ssr"), Some(&"pure_ssr"));
    assert_eq!(by_id.get("r_sub"), Some(&"tr_with_subrepeat"));
    assert_eq!(by_id.get("r_sub_ssr"), Some(&"tr_with_subrepeat_with_ssr"));
}
