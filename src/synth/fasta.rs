//! FASTA writer for the synth pipeline.
//!
//! One record per call, 80-char line width.

use anyhow::{Context, Result};
use std::io::Write;
use std::path::Path;

pub fn write(prefix: &Path, array_id: &str, sequence: &[u8]) -> Result<()> {
    let path = with_ext(prefix, "fa");
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut f = std::fs::File::create(&path).with_context(|| format!("creating {:?}", path))?;
    writeln!(f, ">{}", array_id)?;
    for chunk in sequence.chunks(80) {
        f.write_all(chunk)?;
        f.write_all(b"\n")?;
    }
    Ok(())
}

pub(crate) fn with_ext(prefix: &Path, ext: &str) -> std::path::PathBuf {
    let mut s = prefix.as_os_str().to_owned();
    s.push(".");
    s.push(ext);
    std::path::PathBuf::from(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_record_and_wraps_at_80() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("foo");
        let seq: Vec<u8> = (0..200).map(|i| b"ACGT"[i % 4]).collect();
        write(&prefix, "arr1", &seq).unwrap();
        let text = std::fs::read_to_string(with_ext(&prefix, "fa")).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], ">arr1");
        assert_eq!(lines[1].len(), 80);
        assert_eq!(lines[2].len(), 80);
        assert_eq!(lines[3].len(), 40);
        // No trailing record.
        assert_eq!(lines.len(), 4);
    }
}
