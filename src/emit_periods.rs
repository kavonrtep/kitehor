//! Bridge from `kite-periodicity` output to the v2 detector's
//! `periods.tsv` schema.
//!
//! Score mapping (settled 2026-05-16 in the integration discussion):
//!
//! | rule verdict       | rows emitted                                  |
//! |--------------------|-----------------------------------------------|
//! | `Hor{founder,tile}`| founder @ 0.95 (kite_founder), tile @ 0.90    |
//! |                    |   (kite_tile) if distinct, plus other top-3   |
//! |                    |   peaks @ 0.60 (kite_secondary)               |
//! | `Tandem{monomer}`  | monomer @ 0.95 (kite_founder), plus other     |
//! |                    |   top-3 peaks @ 0.60 (kite_secondary)         |
//! | `Unresolved`       | top-3 peaks @ 0.50 / 0.40 / 0.30 (kite_peak)  |
//! | `NoSignal`         | no rows                                       |
//! | no classifier      | top-3 peaks @ 0.60 (kite_peak)                |
//!
//! Scores are chosen relative to the detector's
//! `DetectorConfig::strong_period_score` (default 0.85): values ≥
//! that floor can fire HOR-rescue and "multi_block_via_strong" paths;
//! below-floor scores act as hints to the canonical column-IC test
//! only. Source labels are documentation-only — the v2 detector
//! gates on `period_score`, not `source`.

use crate::kite::KiteResult;
use crate::rule::RuleVerdict;
use anyhow::{Context, Result};
use std::collections::HashSet;
use std::io::Write;
use std::path::Path;

/// One row of the v2 detector's `periods.tsv` schema.
#[derive(Debug, Clone, PartialEq)]
pub struct PeriodsRow {
    pub array_id: String,
    pub period_bp: usize,
    pub period_score: f64,
    pub source: String,
}

/// Header used by the v2 detector loader (`detect::io::load_periods`).
pub const PERIODS_HEADER: &str = "array_id\tperiod_bp\tperiod_score\tsource";

/// How many non-founder/tile peaks to admit as secondaries / hints.
const MAX_SECONDARIES: usize = 3;

/// Score floor under which everything stays a hint (below
/// detector's default `strong_period_score = 0.85`).
const HINT_SCORE_DESCENDING: [f64; 3] = [0.50, 0.40, 0.30];

/// Score awarded to the rule classifier's founder.
const FOUNDER_SCORE: f64 = 0.95;
/// Score awarded to the rule classifier's tile when it's distinct
/// from the founder.
const TILE_SCORE: f64 = 0.90;
/// Score awarded to other top-N kite peaks once a founder/tile has
/// been identified — visible to the detector but below the rescue
/// threshold.
const SECONDARY_SCORE: f64 = 0.60;

/// Map a `(KiteResult, optional verdict)` pair to v2 detector
/// periods rows.
pub fn build_rows(kr: &KiteResult, verdict: Option<&RuleVerdict>) -> Vec<PeriodsRow> {
    let array_id = kr.array_id.clone();
    let mut rows: Vec<PeriodsRow> = Vec::new();
    let mut used: HashSet<usize> = HashSet::new();

    match verdict {
        Some(RuleVerdict::Hor { founder, tile, .. }) => {
            rows.push(PeriodsRow {
                array_id: array_id.clone(),
                period_bp: *founder,
                period_score: FOUNDER_SCORE,
                source: "kite_founder".into(),
            });
            used.insert(*founder);
            if *tile != *founder {
                rows.push(PeriodsRow {
                    array_id: array_id.clone(),
                    period_bp: *tile,
                    period_score: TILE_SCORE,
                    source: "kite_tile".into(),
                });
                used.insert(*tile);
            }
            append_secondaries(&mut rows, &mut used, kr, &array_id);
        }
        Some(RuleVerdict::Tandem { monomer_bp }) => {
            rows.push(PeriodsRow {
                array_id: array_id.clone(),
                period_bp: *monomer_bp,
                period_score: FOUNDER_SCORE,
                source: "kite_founder".into(),
            });
            used.insert(*monomer_bp);
            append_secondaries(&mut rows, &mut used, kr, &array_id);
        }
        Some(RuleVerdict::Unresolved) => {
            append_hint_peaks(&mut rows, kr, &array_id);
        }
        Some(RuleVerdict::NoSignal) => {
            // No rows. Detector with --allow-missing-periods → ambiguous.
        }
        None => {
            // No classifier ran: emit raw kite peaks as hints at a
            // single below-floor score (no way to tell founder from
            // harmonic without the rule).
            for p in kr.peaks.iter().take(MAX_SECONDARIES) {
                rows.push(PeriodsRow {
                    array_id: array_id.clone(),
                    period_bp: p.period,
                    period_score: SECONDARY_SCORE,
                    source: "kite_peak".into(),
                });
            }
        }
    }
    rows
}

fn append_secondaries(
    rows: &mut Vec<PeriodsRow>,
    used: &mut HashSet<usize>,
    kr: &KiteResult,
    array_id: &str,
) {
    let mut n = 0usize;
    for p in &kr.peaks {
        if n >= MAX_SECONDARIES {
            break;
        }
        if used.contains(&p.period) {
            continue;
        }
        rows.push(PeriodsRow {
            array_id: array_id.to_string(),
            period_bp: p.period,
            period_score: SECONDARY_SCORE,
            source: "kite_secondary".into(),
        });
        used.insert(p.period);
        n += 1;
    }
}

fn append_hint_peaks(rows: &mut Vec<PeriodsRow>, kr: &KiteResult, array_id: &str) {
    for (i, p) in kr.peaks.iter().take(HINT_SCORE_DESCENDING.len()).enumerate() {
        rows.push(PeriodsRow {
            array_id: array_id.to_string(),
            period_bp: p.period,
            period_score: HINT_SCORE_DESCENDING[i],
            source: "kite_peak".into(),
        });
    }
}

/// Write all rows for many arrays to a single `periods.tsv` file.
pub fn write_tsv(path: &Path, batches: &[Vec<PeriodsRow>]) -> Result<usize> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut f = std::fs::File::create(path)
        .with_context(|| format!("creating {:?}", path))?;
    writeln!(f, "{}", PERIODS_HEADER)?;
    let mut n = 0usize;
    for batch in batches {
        for r in batch {
            writeln!(
                f,
                "{}\t{}\t{:.4}\t{}",
                r.array_id, r.period_bp, r.period_score, r.source
            )?;
            n += 1;
        }
    }
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kite::KitePeak;

    fn mk_peak(period: usize, score: f64) -> KitePeak {
        KitePeak {
            period,
            peak_height: 1.0,
            score,
            score2: 0.0,
            score2_norm: 0.01,
            background: 0.0,
        }
    }

    fn mk_result(id: &str, peaks: Vec<KitePeak>) -> KiteResult {
        KiteResult {
            array_id: id.to_string(),
            length_bp: 10_000,
            peaks,
            profile: None,
            background: None,
        }
    }

    #[test]
    fn hor_emits_founder_tile_and_secondaries() {
        let kr = mk_result(
            "a1",
            vec![
                mk_peak(2052, 1.0), // d1 / tile
                mk_peak(171, 0.5),  // founder
                mk_peak(342, 0.2),  // secondary
                mk_peak(513, 0.1),  // secondary
            ],
        );
        let v = RuleVerdict::Hor {
            founder: 171,
            tile: 2052,
            k: 12,
            share: 0.5,
        };
        let rows = build_rows(&kr, Some(&v));
        assert_eq!(rows.len(), 4);
        assert_eq!(rows[0].period_bp, 171);
        assert!((rows[0].period_score - 0.95).abs() < 1e-9);
        assert_eq!(rows[0].source, "kite_founder");
        assert_eq!(rows[1].period_bp, 2052);
        assert!((rows[1].period_score - 0.90).abs() < 1e-9);
        assert_eq!(rows[1].source, "kite_tile");
        // Secondaries
        for r in &rows[2..] {
            assert!((r.period_score - 0.60).abs() < 1e-9);
            assert_eq!(r.source, "kite_secondary");
        }
    }

    #[test]
    fn tandem_emits_single_high_score() {
        let kr = mk_result(
            "a1",
            vec![mk_peak(178, 1.0), mk_peak(356, 0.1)],
        );
        let v = RuleVerdict::Tandem { monomer_bp: 178 };
        let rows = build_rows(&kr, Some(&v));
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].period_bp, 178);
        assert!((rows[0].period_score - 0.95).abs() < 1e-9);
        assert_eq!(rows[1].period_bp, 356);
        assert!((rows[1].period_score - 0.60).abs() < 1e-9);
    }

    #[test]
    fn unresolved_emits_below_floor_hints() {
        let kr = mk_result(
            "a1",
            vec![
                mk_peak(10, 1.0),
                mk_peak(7, 0.5),
                mk_peak(5, 0.2),
                mk_peak(3, 0.1),
            ],
        );
        let rows = build_rows(&kr, Some(&RuleVerdict::Unresolved));
        assert_eq!(rows.len(), 3);
        assert!((rows[0].period_score - 0.50).abs() < 1e-9);
        assert!((rows[1].period_score - 0.40).abs() < 1e-9);
        assert!((rows[2].period_score - 0.30).abs() < 1e-9);
        for r in &rows {
            assert_eq!(r.source, "kite_peak");
            assert!(r.period_score < 0.85, "hints must stay below strong_period_score");
        }
    }

    #[test]
    fn no_signal_emits_nothing() {
        let kr = mk_result("a1", vec![]);
        let rows = build_rows(&kr, Some(&RuleVerdict::NoSignal));
        assert!(rows.is_empty());
    }

    #[test]
    fn no_classifier_emits_secondaries_only() {
        let kr = mk_result(
            "a1",
            vec![mk_peak(171, 1.0), mk_peak(342, 0.5), mk_peak(513, 0.2)],
        );
        let rows = build_rows(&kr, None);
        assert_eq!(rows.len(), 3);
        for r in &rows {
            assert_eq!(r.source, "kite_peak");
            assert!((r.period_score - 0.60).abs() < 1e-9);
            assert!(r.period_score < 0.85);
        }
    }

    #[test]
    fn founder_equals_tile_dedup() {
        // Degenerate case: rule says HOR with founder=tile (shouldn't
        // happen in practice but worth covering).
        let kr = mk_result("a1", vec![mk_peak(171, 1.0)]);
        let v = RuleVerdict::Hor {
            founder: 171,
            tile: 171,
            k: 1,
            share: 1.0,
        };
        let rows = build_rows(&kr, Some(&v));
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].source, "kite_founder");
    }

    #[test]
    fn secondaries_capped_at_max() {
        // 10 peaks; only 3 should make it through after founder.
        let peaks: Vec<KitePeak> = (0..10).map(|i| mk_peak(200 + i, 1.0 - 0.05 * i as f64)).collect();
        let kr = mk_result("a1", peaks);
        let v = RuleVerdict::Tandem { monomer_bp: 200 };
        let rows = build_rows(&kr, Some(&v));
        // 1 founder + at most 3 secondaries
        assert!(rows.len() <= 4);
        let n_sec = rows.iter().filter(|r| r.source == "kite_secondary").count();
        assert_eq!(n_sec, MAX_SECONDARIES);
    }

    #[test]
    fn write_tsv_uses_v2_schema() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("test.periods.tsv");
        let kr = mk_result("a1", vec![mk_peak(171, 1.0), mk_peak(342, 0.1)]);
        let rows = build_rows(&kr, Some(&RuleVerdict::Tandem { monomer_bp: 171 }));
        let n = write_tsv(&p, &[rows]).unwrap();
        assert_eq!(n, 2);
        let s = std::fs::read_to_string(&p).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[0], PERIODS_HEADER);
        assert!(lines[1].starts_with("a1\t171\t0.9500\tkite_founder"));
    }
}
