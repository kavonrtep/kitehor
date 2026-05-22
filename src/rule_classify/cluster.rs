//! Period-cluster construction. Single-linkage on relative period gap;
//! per-cluster aggregation by `score2_norm`-weighted period mean.

use super::decide::PeakRow;

#[derive(Debug, Clone)]
pub struct Cluster {
    pub rep_period: f64,
    pub total_score: f64,
    pub n_peaks: usize,
    pub min_rank: u32,
    pub periods: Vec<u32>,
}

/// Single-linkage cluster peaks by relative period gap. Cuts wherever
/// `(p_cur - p_prev) / p_cur > tol` after sorting by period ascending.
///
/// **NOTE** (audit-confirmed gotcha from
/// `tools/rule_proto/rule_proto.py:84`): the divisor of the relative
/// gap is **`p_cur`** (the later peak after the sort), *not* the
/// cluster mean. Easy to misread when porting.
pub fn cluster_peaks(rows: &[PeakRow], tol: f64) -> Vec<Cluster> {
    if rows.is_empty() {
        return Vec::new();
    }
    let mut sorted: Vec<&PeakRow> = rows.iter().collect();
    sorted.sort_by(|a, b| a.period.cmp(&b.period).then(a.rank.cmp(&b.rank)));

    // Build single-linkage groups.
    let mut groups: Vec<Vec<usize>> = vec![vec![0]];
    for i in 1..sorted.len() {
        let last = *groups.last().unwrap().last().unwrap();
        let p_prev = sorted[last].period as f64;
        let p_cur = sorted[i].period as f64;
        let denom = p_cur.max(1.0);
        let rel_gap = (p_cur - p_prev) / denom;
        if rel_gap > tol {
            groups.push(vec![i]);
        } else {
            groups.last_mut().unwrap().push(i);
        }
    }

    let mut out: Vec<Cluster> = Vec::with_capacity(groups.len());
    for g in &groups {
        let weights: Vec<f64> = g.iter().map(|&i| sorted[i].score2_norm).collect();
        let periods: Vec<f64> = g.iter().map(|&i| sorted[i].period as f64).collect();
        let w_sum: f64 = weights.iter().sum();
        let rep = if w_sum > 0.0 {
            let num: f64 = weights.iter().zip(periods.iter()).map(|(w, p)| w * p).sum();
            num / w_sum
        } else {
            // pandas mean of integer-valued floats: arithmetic mean
            periods.iter().sum::<f64>() / periods.len() as f64
        };
        let total_score: f64 = weights.iter().sum();
        let n_peaks = g.len();
        let min_rank = g.iter().map(|&i| sorted[i].rank).min().unwrap_or(0);
        let p_ints: Vec<u32> = g.iter().map(|&i| sorted[i].period as u32).collect();
        out.push(Cluster {
            rep_period: rep,
            total_score,
            n_peaks,
            min_rank,
            periods: p_ints,
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(rank: u32, period: usize, s: f64) -> PeakRow {
        PeakRow {
            rank,
            period,
            score2_norm: s,
        }
    }

    #[test]
    fn empty_input_yields_no_clusters() {
        assert!(cluster_peaks(&[], 0.015).is_empty());
    }

    #[test]
    fn single_peak_one_cluster() {
        let cs = cluster_peaks(&[row(1, 100, 0.5)], 0.015);
        assert_eq!(cs.len(), 1);
        assert!((cs[0].rep_period - 100.0).abs() < 1e-9);
        assert!((cs[0].total_score - 0.5).abs() < 1e-9);
        assert_eq!(cs[0].n_peaks, 1);
        assert_eq!(cs[0].min_rank, 1);
    }

    #[test]
    fn near_peaks_merge_into_one_cluster() {
        // 300, 301, 302 with tol 0.015 (1.5%). gap(300->301)/301 = 0.0033.
        let rows = vec![row(3, 300, 0.1), row(1, 301, 0.4), row(2, 302, 0.2)];
        let cs = cluster_peaks(&rows, 0.015);
        assert_eq!(cs.len(), 1);
        let c = &cs[0];
        // weighted mean: (0.1*300 + 0.4*301 + 0.2*302) / 0.7
        let want = (0.1 * 300.0 + 0.4 * 301.0 + 0.2 * 302.0) / 0.7;
        assert!((c.rep_period - want).abs() < 1e-6);
        assert_eq!(c.min_rank, 1);
    }

    #[test]
    fn far_peaks_split() {
        // 300 and 900 are 200% apart — definitely a split.
        let rows = vec![row(1, 300, 0.4), row(2, 900, 0.3)];
        let cs = cluster_peaks(&rows, 0.015);
        assert_eq!(cs.len(), 2);
    }

    #[test]
    fn zero_weight_falls_back_to_mean() {
        let rows = vec![row(1, 100, 0.0), row(2, 102, 0.0)];
        let cs = cluster_peaks(&rows, 0.05);
        assert_eq!(cs.len(), 1);
        assert!((cs[0].rep_period - 101.0).abs() < 1e-9);
    }
}
