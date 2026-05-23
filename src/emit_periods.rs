//! Bridge from `kite-periodicity` output to the v2 detector's
//! `periods.tsv` schema.
//!
//! Score mapping (settled 2026-05-16 in the integration discussion;
//! tightened per `docs/reviews/kite_emit_periods_integration_review_2026-05-16.md`
//! review finding #1):
//!
//! | rule verdict       | rows emitted                                  |
//! |--------------------|-----------------------------------------------|
//! | `Hor{founder,tile}`| founder @ 0.95 (`kite_founder`); tile @ 0.90  |
//! |                    |   (`kite_tile`) if ≠ founder; any of the      |
//! |                    |   **top-3 Kite peaks** not already used @ 0.60|
//! |                    |   (`kite_secondary`)                          |
//! | `Tandem{monomer}`  | monomer @ 0.95 (`kite_monomer`); any of the   |
//! |                    |   **top-3 Kite peaks** not already used @ 0.60|
//! |                    |   (`kite_secondary`)                          |
//! | `Unresolved`       | top-3 peaks @ 0.50 / 0.40 / 0.30 (`kite_peak`)|
//! | `NoSignal`         | no rows                                       |
//! | no classifier      | top-3 peaks @ 0.60 (`kite_peak`)              |
//!
//! **Top-3 contract**: the emitter never looks past `kr.peaks[0..3]`.
//! Founder/tile peaks already in the top-3 don't double up;
//! founder/tile peaks outside the top-3 (rare; the rule classifier
//! requires founder in top-N where N is configurable upstream) are
//! kept as the high-score row and don't earn an additional
//! secondary slot.
//!
//! Scores are chosen relative to the detector's
//! `DetectorConfig::strong_period_score` (default 0.85): values ≥
//! that floor can fire HOR-rescue and "multi_block_via_strong" paths;
//! below-floor scores act as hints to the canonical column-IC test
//! only. Source labels are documentation-only — the v2 detector
//! gates on `period_score`, not `source`.

use crate::kite::KiteResult;
use crate::rule_classify::LegacyVerdict as RuleVerdict;
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

/// Size of the Kite-peak window the emitter considers. Review-#1:
/// secondaries are drawn only from the top-3, NOT "the next three
/// peaks after excluding founder/tile". Hint-only paths also cap at 3.
const TOP_N_KITE_PEAKS: usize = 3;

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
            // Review-#5: rename for clarity. A tandem monomer is not a
            // HOR founder. Detector ignores `source`, so this is a
            // user-facing relabel only.
            rows.push(PeriodsRow {
                array_id: array_id.clone(),
                period_bp: *monomer_bp,
                period_score: FOUNDER_SCORE,
                source: "kite_monomer".into(),
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
            for p in kr.peaks.iter().take(TOP_N_KITE_PEAKS) {
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
    // Review-#1: only walk Kite's top-3 peaks. Anything past rank 3
    // is not part of the documented contract and risks letting
    // harmonic / noisy lower-rank peaks reach the detector.
    for p in kr.peaks.iter().take(TOP_N_KITE_PEAKS) {
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
    }
}

fn append_hint_peaks(rows: &mut Vec<PeriodsRow>, kr: &KiteResult, array_id: &str) {
    let n = TOP_N_KITE_PEAKS.min(HINT_SCORE_DESCENDING.len());
    for (i, p) in kr.peaks.iter().take(n).enumerate() {
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
    let mut f = std::fs::File::create(path).with_context(|| format!("creating {:?}", path))?;
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
            profile: Vec::new(),
            background: Vec::new(),
        }
    }

    #[test]
    fn hor_emits_founder_tile_and_secondaries() {
        // Top-3 are [2052, 171, 342]; founder=171, tile=2052 → only
        // 342 qualifies as a secondary. Rank-4 (513) is excluded by
        // the top-3 cap (Review-#1).
        let kr = mk_result(
            "a1",
            vec![
                mk_peak(2052, 1.0), // d1 / tile
                mk_peak(171, 0.5),  // founder
                mk_peak(342, 0.2),  // secondary
                mk_peak(513, 0.1),  // rank 4 — excluded
            ],
        );
        let v = RuleVerdict::Hor {
            founder: 171,
            tile: 2052,
            k: 12,
            share: 0.5,
        };
        let rows = build_rows(&kr, Some(&v));
        assert_eq!(rows.len(), 3, "founder + tile + 1 top-3-eligible secondary");
        assert_eq!(rows[0].period_bp, 171);
        assert!((rows[0].period_score - 0.95).abs() < 1e-9);
        assert_eq!(rows[0].source, "kite_founder");
        assert_eq!(rows[1].period_bp, 2052);
        assert!((rows[1].period_score - 0.90).abs() < 1e-9);
        assert_eq!(rows[1].source, "kite_tile");
        // The single secondary.
        assert_eq!(rows[2].period_bp, 342);
        assert!((rows[2].period_score - 0.60).abs() < 1e-9);
        assert_eq!(rows[2].source, "kite_secondary");
    }

    #[test]
    fn tandem_emits_single_high_score() {
        let kr = mk_result("a1", vec![mk_peak(178, 1.0), mk_peak(356, 0.1)]);
        let v = RuleVerdict::Tandem { monomer_bp: 178 };
        let rows = build_rows(&kr, Some(&v));
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].period_bp, 178);
        assert!((rows[0].period_score - 0.95).abs() < 1e-9);
        // Review-#5: tandem monomer is labeled `kite_monomer`,
        // not `kite_founder`.
        assert_eq!(rows[0].source, "kite_monomer");
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
            assert!(
                r.period_score < 0.85,
                "hints must stay below strong_period_score"
            );
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
    fn secondaries_capped_at_top_3() {
        // 10 peaks; monomer is rank 1 (200). Top-3 are [200, 201, 202];
        // monomer takes 200 → secondaries = {201, 202}.
        let peaks: Vec<KitePeak> = (0..10)
            .map(|i| mk_peak(200 + i, 1.0 - 0.05 * i as f64))
            .collect();
        let kr = mk_result("a1", peaks);
        let v = RuleVerdict::Tandem { monomer_bp: 200 };
        let rows = build_rows(&kr, Some(&v));
        // 1 monomer + 2 secondaries (the two remaining top-3 peaks).
        assert_eq!(rows.len(), 3);
        let n_sec = rows.iter().filter(|r| r.source == "kite_secondary").count();
        assert_eq!(n_sec, 2);
        let sec_periods: Vec<usize> = rows
            .iter()
            .filter(|r| r.source == "kite_secondary")
            .map(|r| r.period_bp)
            .collect();
        assert_eq!(sec_periods, vec![201, 202]);
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
        // Review-#5: tandem high-score row labels source `kite_monomer`.
        assert!(lines[1].starts_with("a1\t171\t0.9500\tkite_monomer"));
    }

    // Review-#1: secondaries must not include peaks past rank 3.
    #[test]
    fn secondaries_never_look_past_top_3() {
        // 5 peaks; founder at rank 5 (intentionally outside top-3 —
        // wouldn't actually happen via the rule classifier but
        // exercises the cap). Secondaries are drawn from top-3 only,
        // so we should see at most founder + 3 secondaries.
        let kr = mk_result(
            "a1",
            vec![
                mk_peak(100, 1.0), // rank 1
                mk_peak(200, 0.8), // rank 2
                mk_peak(300, 0.6), // rank 3
                mk_peak(400, 0.4), // rank 4 — never emitted
                mk_peak(500, 0.2), // rank 5 — founder for this test
            ],
        );
        let v = RuleVerdict::Hor {
            founder: 500,
            tile: 1000, // tile not in peaks list
            k: 2,
            share: 0.5,
        };
        let rows = build_rows(&kr, Some(&v));
        let periods: Vec<usize> = rows.iter().map(|r| r.period_bp).collect();
        assert!(periods.contains(&500));
        assert!(periods.contains(&1000));
        // None of the rank-4 or rank-5 (= founder, allowed) periods
        // should appear as a secondary. Rank 4 (400) must be absent.
        assert!(!periods.contains(&400), "rank-4 peak leaked past top-3 cap");
        // Top-3 peaks 100/200/300 each should appear as secondaries.
        for w in [100, 200, 300] {
            assert!(periods.contains(&w), "top-3 peak {w} should be a secondary");
        }
    }

    #[test]
    fn secondaries_skip_founder_and_tile_in_top_3() {
        // founder + tile both fall within Kite's top-3 → they should
        // not double as secondaries. Only the remaining top-3 peak
        // should fire as a secondary.
        let kr = mk_result(
            "a1",
            vec![
                mk_peak(2052, 1.0), // tile, rank 1
                mk_peak(171, 0.8),  // founder, rank 2
                mk_peak(342, 0.6),  // the one valid secondary
                mk_peak(513, 0.4),  // rank 4 — must NOT be emitted
            ],
        );
        let v = RuleVerdict::Hor {
            founder: 171,
            tile: 2052,
            k: 12,
            share: 0.5,
        };
        let rows = build_rows(&kr, Some(&v));
        let secondaries: Vec<&PeriodsRow> = rows
            .iter()
            .filter(|r| r.source == "kite_secondary")
            .collect();
        assert_eq!(secondaries.len(), 1);
        assert_eq!(secondaries[0].period_bp, 342);
    }
}
