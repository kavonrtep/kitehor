//! Peaks-TSV reader/writer for `rescore`.
//!
//! We parse the header to map column names → indices (so a future column
//! reorder upstream doesn't silently corrupt rescore output), then keep
//! each row's original text verbatim so the writer can append columns
//! without reformatting the upstream floats — preserving byte equality on
//! the unchanged cells.

use anyhow::{anyhow, Context, Result};
use std::fs;
use std::path::Path;

/// Indices (0-based) of the columns we need to read out of every row.
/// Other columns are preserved unchanged in the row's verbatim text.
#[derive(Debug, Clone, Copy)]
struct ColumnIdx {
    case_id: usize,
    rank: usize,
    period: usize,
}

fn parse_header(header: &str) -> Result<ColumnIdx> {
    let cols: Vec<&str> = header.split('\t').collect();
    let find = |name: &str| -> Result<usize> {
        cols.iter()
            .position(|c| *c == name)
            .ok_or_else(|| anyhow!("peaks header missing required column: {}", name))
    };
    Ok(ColumnIdx {
        case_id: find("case_id")?,
        rank: find("rank")?,
        period: find("period")?,
    })
}

/// One row from the input peaks TSV.
#[derive(Debug, Clone)]
pub struct RawRow {
    /// Original line, unchanged (tab-separated, no trailing newline).
    pub line: String,
    pub case_id: String,
    pub rank: usize,
    pub period: usize,
}

/// Parsed peaks file: header line + ordered rows.
#[derive(Debug, Clone)]
pub struct LoadedPeaks {
    pub header: String,
    pub rows: Vec<RawRow>,
}

pub fn load_peaks(path: &Path) -> Result<LoadedPeaks> {
    let content =
        fs::read_to_string(path).with_context(|| format!("reading peaks file {:?}", path))?;
    let mut iter = content.lines();
    let header = iter
        .next()
        .ok_or_else(|| anyhow!("peaks file {:?} is empty", path))?
        .to_string();
    let idx = parse_header(&header)?;

    let mut rows = Vec::new();
    for (i, line) in iter.enumerate() {
        if line.is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        let case_id = cells
            .get(idx.case_id)
            .ok_or_else(|| anyhow!("row {} short of case_id", i + 2))?
            .to_string();
        let rank: usize = cells
            .get(idx.rank)
            .ok_or_else(|| anyhow!("row {} short of rank", i + 2))?
            .parse()
            .with_context(|| format!("row {}: parsing rank", i + 2))?;
        let period: usize = cells
            .get(idx.period)
            .ok_or_else(|| anyhow!("row {} short of period", i + 2))?
            .parse()
            .with_context(|| format!("row {}: parsing period", i + 2))?;
        rows.push(RawRow {
            line: line.to_string(),
            case_id,
            rank,
            period,
        });
    }
    Ok(LoadedPeaks { header, rows })
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn parses_kite_header() {
        let h = "case_id\tarray_length\trank\tperiod\tpeak_height\tscore\tscore2\tscore2_norm\tbackground";
        let idx = parse_header(h).unwrap();
        assert_eq!(idx.case_id, 0);
        assert_eq!(idx.rank, 2);
        assert_eq!(idx.period, 3);
    }

    #[test]
    fn missing_column_errors() {
        assert!(parse_header("case_id\trank").is_err());
    }

    #[test]
    fn load_peaks_roundtrip() {
        let mut f = NamedTempFile::new().unwrap();
        use std::io::Write;
        writeln!(
            f,
            "case_id\tarray_length\trank\tperiod\tpeak_height\tscore\tscore2\tscore2_norm\tbackground"
        )
        .unwrap();
        writeln!(f, "rec1\t1000\t1\t200\t100.0\t0.5\t1.0\t0.9\t10.0").unwrap();
        writeln!(f, "rec1\t1000\t2\t300\t80.0\t0.4\t0.8\t0.7\t10.0").unwrap();
        let loaded = load_peaks(f.path()).unwrap();
        assert_eq!(loaded.rows.len(), 2);
        assert_eq!(loaded.rows[0].case_id, "rec1");
        assert_eq!(loaded.rows[0].rank, 1);
        assert_eq!(loaded.rows[0].period, 200);
        assert_eq!(loaded.rows[1].period, 300);
        // Verbatim line preservation.
        assert_eq!(
            loaded.rows[0].line,
            "rec1\t1000\t1\t200\t100.0\t0.5\t1.0\t0.9\t10.0"
        );
    }
}
