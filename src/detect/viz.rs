//! Visualisation + matrix export (`detect_impl_plan.md §8`, A7).
//!
//! Per array, into `<viz_dir>/<array_id>/`:
//!
//! - `raster_w{w}.tsv`            — numeric wrapped matrix (0=A, 1=C, 2=G, 3=T, 4=N)
//! - `raster_w{w}.png`            — same matrix, A/C/G/T colour-coded (feature `viz`)
//! - `column_ic_w{w}.tsv`         — per-column IC at width w
//! - `rk_w{w}.tsv`                — R(k) curve for k = 1..K
//! - `shift_w{w}.tsv`             — best_shift(r) signal
//! - `column_edge_rate_w{w}.tsv`  — diff_y per column
//! - `edge_matrix_w{w}.tsv`       — per-row diff_y bits (one row per (r, r+1) pair)
//!
//! Review-2026-05-16 #9 — flag semantics, made exact:
//!
//! - `--viz-dir` alone (no granular flags): writes every cheap TSV
//!   (`column_ic`, `column_edge_rate`, `rk`, `shift`). This is the
//!   documented back-compat default.
//! - With any granular `--export-*` flag set, only the flagged
//!   artefacts are written: each flag gates exactly one output.
//!     - `--export-ic`     → `column_ic_w{w}.tsv`
//!     - `--export-shift`  → `shift_w{w}.tsv` + `rk_w{w}.tsv`
//!     - `--export-edges`  → `column_edge_rate_w{w}.tsv` + `edge_matrix_w{w}.tsv`
//!     - `--export-raster` → `raster_w{w}.tsv` + `raster_w{w}.png`
//!
//! PNG requires the `viz` Cargo feature (default-on). A no-viz build
//! accepts the flags at the CLI surface but returns a clear runtime
//! error if a PNG path actually fires (per A7).

use anyhow::{Context, Result};
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Configuration knobs supplied by the CLI.
#[derive(Debug, Clone, Default)]
pub struct VizFlags {
    pub viz_dir: Option<PathBuf>,
    pub export_raster: bool,
    pub export_shift: bool,
    pub export_edges: bool,
    pub export_ic: bool,
}

impl VizFlags {
    pub fn is_active(&self) -> bool {
        self.viz_dir.is_some()
            || self.export_raster
            || self.export_shift
            || self.export_edges
            || self.export_ic
    }

    /// The output directory; when `viz_dir` is unset but `--export-*`
    /// flags are, fall back to `<cwd>/synth_viz`.
    pub fn dir(&self) -> Option<PathBuf> {
        self.viz_dir.clone()
    }

    /// True iff any granular flag was explicitly set. When false,
    /// `export()` falls back to the "all cheap TSVs" default.
    /// Review-2026-05-16 #9.
    pub fn granular_set(&self) -> bool {
        self.export_raster || self.export_shift || self.export_edges || self.export_ic
    }
}

/// Per-array bundle of inputs the viz layer needs. Builders inside
/// `detect::mod::run_one` populate as much as is available.
#[derive(Debug, Clone)]
pub struct VizBundle<'a> {
    pub array_id: &'a str,
    pub width_bp: usize,
    pub seq: &'a [u8],
    pub n_rows: usize,
    pub column_ic: Option<&'a [f64]>,
    pub column_edge_rate: Option<&'a [f64]>,
    pub r_k: Option<&'a [f64]>,
    pub best_shift: Option<&'a [i32]>,
}

pub fn export(flags: &VizFlags, bundle: &VizBundle<'_>) -> Result<()> {
    let Some(root) = flags.dir() else {
        return Ok(());
    };
    let dir = root.join(bundle.array_id);
    std::fs::create_dir_all(&dir)
        .with_context(|| format!("creating viz dir {:?}", dir))?;

    // Review-2026-05-16 #9: each granular flag gates exactly one
    // (or one bundle of) artefact. `--viz-dir` alone (no granular
    // flags) → emit every cheap TSV (back-compat default).
    let granular = flags.granular_set();
    let emit_ic = !granular || flags.export_ic;
    let emit_shift = !granular || flags.export_shift;
    let emit_edges = !granular || flags.export_edges;

    if emit_ic {
        if let Some(ic) = bundle.column_ic {
            write_one_d_tsv(&dir.join(format!("column_ic_w{}.tsv", bundle.width_bp)), ic, "ic")?;
        }
    }
    if emit_edges {
        if let Some(rate) = bundle.column_edge_rate {
            write_one_d_tsv(
                &dir.join(format!("column_edge_rate_w{}.tsv", bundle.width_bp)),
                rate,
                "edge_rate",
            )?;
        }
        // `--export-edges` (explicit) also writes the per-row diff_y
        // matrix the reviewer flagged as missing. Skipped in the
        // "viz_dir alone" default so the back-compat output set
        // stays cheap.
        if flags.export_edges && bundle.n_rows >= 2 {
            write_edge_matrix_tsv(
                &dir.join(format!("edge_matrix_w{}.tsv", bundle.width_bp)),
                bundle.seq,
                bundle.width_bp,
                bundle.n_rows,
            )?;
        }
    }
    if emit_shift {
        if let Some(r) = bundle.r_k {
            write_one_d_tsv(&dir.join(format!("rk_w{}.tsv", bundle.width_bp)), r, "r_k")?;
        }
        if let Some(s) = bundle.best_shift {
            write_shift_tsv(
                &dir.join(format!("shift_w{}.tsv", bundle.width_bp)),
                s,
            )?;
        }
    }

    if flags.export_raster {
        write_raster_tsv(
            &dir.join(format!("raster_w{}.tsv", bundle.width_bp)),
            bundle.seq,
            bundle.width_bp,
            bundle.n_rows,
        )?;
        write_raster_png(
            &dir.join(format!("raster_w{}.png", bundle.width_bp)),
            bundle.seq,
            bundle.width_bp,
            bundle.n_rows,
        )?;
    }
    Ok(())
}

/// Per-row diff_y matrix: one row per (r, r+1) pair, one column per
/// position c. `1` = bases differ, `0` = bases match. Empty if
/// `n_rows < 2`.
fn write_edge_matrix_tsv(path: &Path, seq: &[u8], width: usize, n_rows: usize) -> Result<()> {
    let mut f = File::create(path).with_context(|| format!("creating {:?}", path))?;
    write!(f, "# array width_bp={} n_pairs={} schema_version=1\npair", width, n_rows - 1)?;
    for c in 0..width {
        write!(f, "\tcol_{}", c)?;
    }
    writeln!(f)?;
    for r in 0..n_rows - 1 {
        write!(f, "{}", r)?;
        let base = r * width;
        let next = (r + 1) * width;
        for c in 0..width {
            let bit = (seq[base + c] != seq[next + c]) as u8;
            write!(f, "\t{bit}")?;
        }
        writeln!(f)?;
    }
    Ok(())
}

fn write_one_d_tsv(path: &Path, xs: &[f64], col_name: &str) -> Result<()> {
    let mut f = File::create(path).with_context(|| format!("creating {:?}", path))?;
    writeln!(f, "index\t{col_name}")?;
    for (i, v) in xs.iter().enumerate() {
        writeln!(f, "{i}\t{:.6}", v)?;
    }
    Ok(())
}

fn write_shift_tsv(path: &Path, xs: &[i32]) -> Result<()> {
    let mut f = File::create(path).with_context(|| format!("creating {:?}", path))?;
    writeln!(f, "row\tbest_shift_bp")?;
    for (i, v) in xs.iter().enumerate() {
        writeln!(f, "{i}\t{v}")?;
    }
    Ok(())
}

/// Numeric wrapped matrix (A=0, C=1, G=2, T=3, N=4). One row per
/// wrap row; tab-separated columns. Useful for downstream Python
/// analysis without parsing FASTA.
pub fn write_raster_tsv(path: &Path, seq: &[u8], width: usize, n_rows: usize) -> Result<()> {
    let mut f = File::create(path).with_context(|| format!("creating {:?}", path))?;
    write!(f, "# array width_bp={} n_rows={} schema_version=1\nrow", width, n_rows)?;
    for c in 0..width {
        write!(f, "\tcol_{}", c)?;
    }
    writeln!(f)?;
    for r in 0..n_rows {
        write!(f, "{}", r)?;
        for c in 0..width {
            let code = base_code(seq[r * width + c]);
            write!(f, "\t{code}")?;
        }
        writeln!(f)?;
    }
    Ok(())
}

#[inline]
fn base_code(b: u8) -> u8 {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

// ---------------- PNG export ----------------

#[cfg(feature = "viz")]
fn write_raster_png(path: &Path, seq: &[u8], width: usize, n_rows: usize) -> Result<()> {
    // Per-base colour palette (matplotlib defaults).
    let colour = |b: u8| -> [u8; 3] {
        match b {
            b'A' => [0x2c, 0xa0, 0x2c],
            b'C' => [0x1f, 0x77, 0xb4],
            b'G' => [0xd6, 0x27, 0x28],
            b'T' => [0x94, 0x67, 0xbd],
            _ => [0x99, 0x99, 0x99],
        }
    };

    // Vertical downsample if more than 4096 rows so PNGs stay sane.
    const ROW_CAP: usize = 4096;
    let stride = (n_rows + ROW_CAP - 1) / ROW_CAP;
    let img_rows = (n_rows + stride - 1) / stride;
    let img_w = u32::try_from(width).context("width too large for PNG")?;
    let img_h = u32::try_from(img_rows).context("img_rows too large for PNG")?;
    let mut buf: Vec<u8> = Vec::with_capacity(width * img_rows * 3);

    for r in 0..img_rows {
        // Pick the first row of each downsample block (cheap; users
        // who care about averaging can use the TSV).
        let src_r = r * stride;
        for c in 0..width {
            let rgb = colour(seq[src_r * width + c]);
            buf.extend_from_slice(&rgb);
        }
    }
    let img = image::RgbImage::from_raw(img_w, img_h, buf)
        .ok_or_else(|| anyhow::anyhow!("PNG buffer size mismatch"))?;
    img.save(path)
        .with_context(|| format!("saving PNG {:?}", path))?;
    Ok(())
}

#[cfg(not(feature = "viz"))]
fn write_raster_png(_path: &Path, _seq: &[u8], _width: usize, _n_rows: usize) -> Result<()> {
    anyhow::bail!(
        "PNG visualisation support was not compiled in (`viz` Cargo feature disabled). \
         Re-build with `--features viz` or drop --export-raster. TSV diagnostics still \
         work without re-building."
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tsv_is_written_under_array_id_subdir() {
        let dir = tempfile::tempdir().unwrap();
        let flags = VizFlags {
            viz_dir: Some(dir.path().to_path_buf()),
            ..Default::default()
        };
        let ic = [1.0, 2.0, 0.5];
        let bundle = VizBundle {
            array_id: "arr1",
            width_bp: 171,
            seq: b"ACGT",
            n_rows: 1,
            column_ic: Some(&ic),
            column_edge_rate: None,
            r_k: None,
            best_shift: None,
        };
        export(&flags, &bundle).unwrap();
        let p = dir.path().join("arr1").join("column_ic_w171.tsv");
        let s = std::fs::read_to_string(&p).unwrap();
        assert!(s.starts_with("index\tic"));
        assert!(s.contains("0\t1.000000"));
        assert!(s.contains("2\t0.500000"));
    }

    #[test]
    fn raster_tsv_encodes_bases_to_integers() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("r.tsv");
        write_raster_tsv(&p, b"ACGTACGT", 4, 2).unwrap();
        let s = std::fs::read_to_string(&p).unwrap();
        // header + 2 data rows.
        assert_eq!(s.lines().count(), 4);
        let last = s.lines().last().unwrap();
        // row 1: ACGT → 0 1 2 3
        assert!(last.starts_with("1\t0\t1\t2\t3"));
    }

    #[test]
    fn export_off_when_no_flags() {
        let flags = VizFlags::default();
        let ic = [1.0];
        let bundle = VizBundle {
            array_id: "arr",
            width_bp: 100,
            seq: b"A",
            n_rows: 0,
            column_ic: Some(&ic),
            column_edge_rate: None,
            r_k: None,
            best_shift: None,
        };
        // No viz_dir means no I/O; just succeeds.
        export(&flags, &bundle).unwrap();
    }

    // Review-2026-05-16 #9: granular flags must gate their output.
    #[test]
    fn granular_export_only_writes_flagged_artifacts() {
        let dir = tempfile::tempdir().unwrap();
        let flags = VizFlags {
            viz_dir: Some(dir.path().to_path_buf()),
            export_ic: true,
            ..Default::default()
        };
        let ic = [1.0, 2.0, 0.5];
        let rate = [0.1, 0.2, 0.3];
        let rk = [0.0, 0.8];
        let shift = [0i32, 1];
        let bundle = VizBundle {
            array_id: "arr1",
            width_bp: 8,
            seq: b"ACGTACGT",
            n_rows: 1,
            column_ic: Some(&ic),
            column_edge_rate: Some(&rate),
            r_k: Some(&rk),
            best_shift: Some(&shift),
        };
        export(&flags, &bundle).unwrap();
        let base = dir.path().join("arr1");
        assert!(base.join("column_ic_w8.tsv").exists(), "IC was flagged on");
        assert!(
            !base.join("column_edge_rate_w8.tsv").exists(),
            "edge_rate was NOT flagged; should be absent under granular mode"
        );
        assert!(
            !base.join("rk_w8.tsv").exists(),
            "rk was NOT flagged"
        );
        assert!(
            !base.join("shift_w8.tsv").exists(),
            "shift was NOT flagged"
        );
    }

    #[test]
    fn export_edges_writes_edge_matrix() {
        let dir = tempfile::tempdir().unwrap();
        let flags = VizFlags {
            viz_dir: Some(dir.path().to_path_buf()),
            export_edges: true,
            ..Default::default()
        };
        let rate = [0.5, 0.0, 1.0, 0.5];
        let bundle = VizBundle {
            array_id: "arr1",
            width_bp: 4,
            seq: b"ACGTAGGT", // 2 rows: ACGT vs AGGT → diffs at col 1,2
            n_rows: 2,
            column_ic: None,
            column_edge_rate: Some(&rate),
            r_k: None,
            best_shift: None,
        };
        export(&flags, &bundle).unwrap();
        let p = dir.path().join("arr1").join("edge_matrix_w4.tsv");
        let s = std::fs::read_to_string(&p).unwrap();
        // Header comment + column header + 1 pair row.
        assert_eq!(s.lines().count(), 3);
        let last = s.lines().last().unwrap();
        // pair 0: ACGT vs AGGT → 0,1,1,0 at col 0..3
        assert_eq!(last, "0\t0\t1\t1\t0");
    }

    #[cfg(feature = "viz")]
    #[test]
    fn png_writes_a_valid_image() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("r.png");
        write_raster_png(&p, b"ACGTACGTACGTACGT", 4, 4).unwrap();
        let bytes = std::fs::read(&p).unwrap();
        // Standard PNG signature.
        assert_eq!(&bytes[..8], b"\x89PNG\r\n\x1a\n");
    }
}
