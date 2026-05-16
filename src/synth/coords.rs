//! Coord map: logical (block, copy, slot) → realised bp range.
//!
//! Built incrementally by `blocks::expand`. Later stages (`wobble`,
//! `events`, `noise`) update it as they insert/delete bases.
//!
//! MVP representation: a Vec of entries plus a linear-scan `find`.
//! For typical inputs (a few Mb, <100 events) this is well under 1 s;
//! switch to a hash or interval tree only if profiling demands.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CoordEntry {
    pub block_idx: usize,
    /// 1-indexed within the targeted block (matches the YAML schema).
    pub copy_idx: usize,
    /// 1-indexed slot within an HOR copy. Always 1 for SIMPLE_TR.
    pub slot_idx: usize,
    pub realised_start_bp: usize,
    pub realised_len_bp: usize,
}

impl CoordEntry {
    pub fn end_bp(&self) -> usize {
        self.realised_start_bp + self.realised_len_bp
    }
}

#[derive(Debug, Clone, Default)]
pub struct CoordMap {
    pub entries: Vec<CoordEntry>,
}

/// Apply the same indel rule that `CoordMap::apply_indels` uses to a
/// single `(start, len)` span. Returns the updated `(start, len)`.
///
/// Shared with `wobble`, `noise`, and `events` so every stage updates
/// `filler_spans` consistently with `coord_map`. See F3 in the synth
/// implementation review.
pub fn apply_indels_to_span(start: usize, len: usize, indels: &[(usize, i32)]) -> (usize, usize) {
    let s = start;
    let e = s + len;
    let mut shift: i64 = 0;
    let mut len_delta: i64 = 0;
    for &(pos, delta) in indels {
        if pos < s {
            shift += delta as i64;
        } else if pos < e {
            len_delta += delta as i64;
        }
    }
    let new_start = (s as i64 + shift) as usize;
    let new_len = ((len as i64) + len_delta).max(0) as usize;
    (new_start, new_len)
}

/// Shift a `(start, len)` span by `delta` if `start >= pos`. Mirrors
/// `CoordMap::shift_after`.
pub fn shift_span_after(start: usize, len: usize, pos: usize, delta: i64) -> (usize, usize) {
    if start >= pos {
        ((start as i64 + delta) as usize, len)
    } else {
        (start, len)
    }
}

impl CoordMap {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn push(&mut self, e: CoordEntry) {
        self.entries.push(e);
    }
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    /// 1-indexed lookup. Returns the entry matching `(block, copy, slot)`.
    pub fn find(&self, block: usize, copy: usize, slot: usize) -> Option<&CoordEntry> {
        self.entries
            .iter()
            .find(|e| e.block_idx == block && e.copy_idx == copy && e.slot_idx == slot)
    }
    /// All entries in a contiguous copy range of one block, sorted by
    /// `(copy_idx, slot_idx)` ascending.
    pub fn range(&self, block: usize, first_copy: usize, last_copy: usize) -> Vec<&CoordEntry> {
        let mut v: Vec<&CoordEntry> = self
            .entries
            .iter()
            .filter(|e| e.block_idx == block && e.copy_idx >= first_copy && e.copy_idx <= last_copy)
            .collect();
        v.sort_by_key(|e| (e.copy_idx, e.slot_idx));
        v
    }

    /// Shift entries that **start at or after** `pos` by `delta` bp.
    /// Entries strictly before `pos` are untouched (including their length).
    ///
    /// Use this for **uncovered structural insertions** like
    /// `DUPLICATION`: the inserted bytes are owned by an event log entry,
    /// not by any existing coord entry, so the byte immediately after the
    /// insertion (formerly at `pos`) belongs to a downstream entry that
    /// must shift right rather than absorb the insertion.
    ///
    /// This is the policy split flagged in the synth review (F2): the
    /// generic `apply_indels` treats an insertion at `entry.start` as
    /// "inside" the entry (correct for noise/wobble byte-level edits),
    /// but that policy corrupts coord ownership for structural fillers.
    pub fn shift_after(&mut self, pos: usize, delta: i64) {
        for e in &mut self.entries {
            if e.realised_start_bp >= pos {
                e.realised_start_bp = (e.realised_start_bp as i64 + delta) as usize;
            }
        }
    }

    /// Apply a list of `(position, delta)` indels to every entry.
    ///
    /// `delta = +1` is an insertion **at** `position` (the new base
    /// ends up at `position` in the new sequence; original bytes from
    /// `position` onward shift right by 1). `delta = -1` deletes the
    /// byte at `position`.
    ///
    /// Both operations use the **same** classification of `pos`
    /// against each entry — the new/removed byte counts as "part of
    /// the entry" iff `entry.start <= pos < entry.end`:
    ///
    /// - `pos < entry.start`: entry shifts (start += delta).
    /// - `entry.start <= pos < entry.end`: entry length changes (+1
    ///   for insertion, -1 for deletion).
    /// - `pos >= entry.end`: entry unaffected.
    ///
    /// This convention guarantees that when entries fully tile a
    /// contiguous region, the sum of their lengths after `apply_indels`
    /// equals the post-noise length of that region — i.e., every
    /// inserted byte is owned by exactly one entry (the one to its
    /// right at a boundary).
    ///
    /// `indels` must be in **pre-noise** coordinates (positions into
    /// the sequence at the time of the indel, before any subsequent
    /// indel shifts).
    pub fn apply_indels(&mut self, indels: &[(usize, i32)]) {
        for entry in &mut self.entries {
            let s = entry.realised_start_bp;
            let e = s + entry.realised_len_bp;
            let mut shift: i64 = 0;
            let mut len_delta: i64 = 0;
            for &(pos, delta) in indels {
                if pos < s {
                    shift += delta as i64;
                } else if pos < e {
                    len_delta += delta as i64;
                }
            }
            entry.realised_start_bp = (s as i64 + shift) as usize;
            entry.realised_len_bp = (entry.realised_len_bp as i64 + len_delta) as usize;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(b: usize, c: usize, s: usize, start: usize, len: usize) -> CoordEntry {
        CoordEntry {
            block_idx: b,
            copy_idx: c,
            slot_idx: s,
            realised_start_bp: start,
            realised_len_bp: len,
        }
    }

    #[test]
    fn find_returns_matching_entry() {
        let mut m = CoordMap::new();
        m.push(entry(0, 1, 1, 0, 100));
        m.push(entry(0, 1, 2, 100, 100));
        m.push(entry(0, 2, 1, 200, 100));
        let e = m.find(0, 1, 2).unwrap();
        assert_eq!(e.realised_start_bp, 100);
    }

    #[test]
    fn range_filters_and_sorts() {
        let mut m = CoordMap::new();
        m.push(entry(0, 3, 1, 500, 50));
        m.push(entry(0, 1, 1, 0, 50));
        m.push(entry(0, 2, 1, 50, 50));
        m.push(entry(1, 1, 1, 100, 50)); // different block
        let r = m.range(0, 1, 2);
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].copy_idx, 1);
        assert_eq!(r[1].copy_idx, 2);
    }

    #[test]
    fn end_bp_is_start_plus_len() {
        let e = entry(0, 1, 1, 1000, 171);
        assert_eq!(e.end_bp(), 1171);
    }

    #[test]
    fn apply_indels_insertion_before_shifts() {
        let mut m = CoordMap::new();
        m.push(entry(0, 1, 1, 100, 50));
        m.apply_indels(&[(50, 1)]);
        assert_eq!(m.entries[0].realised_start_bp, 101);
        assert_eq!(m.entries[0].realised_len_bp, 50);
    }

    #[test]
    fn apply_indels_insertion_inside_extends() {
        let mut m = CoordMap::new();
        m.push(entry(0, 1, 1, 100, 50));
        m.apply_indels(&[(120, 1)]);
        assert_eq!(m.entries[0].realised_start_bp, 100);
        assert_eq!(m.entries[0].realised_len_bp, 51);
    }

    #[test]
    fn apply_indels_insertion_after_unaffected() {
        let mut m = CoordMap::new();
        m.push(entry(0, 1, 1, 100, 50));
        m.apply_indels(&[(150, 1), (200, 1)]);
        assert_eq!(m.entries[0].realised_start_bp, 100);
        assert_eq!(m.entries[0].realised_len_bp, 50);
    }

    #[test]
    fn apply_indels_deletion_before_shifts_left() {
        let mut m = CoordMap::new();
        m.push(entry(0, 1, 1, 100, 50));
        m.apply_indels(&[(50, -1)]);
        assert_eq!(m.entries[0].realised_start_bp, 99);
        assert_eq!(m.entries[0].realised_len_bp, 50);
    }

    #[test]
    fn apply_indels_deletion_inside_shrinks() {
        let mut m = CoordMap::new();
        m.push(entry(0, 1, 1, 100, 50));
        m.apply_indels(&[(125, -1), (130, -1)]);
        assert_eq!(m.entries[0].realised_start_bp, 100);
        assert_eq!(m.entries[0].realised_len_bp, 48);
    }

    #[test]
    fn apply_indels_at_start_boundary() {
        // Both insertion and deletion at pos == entry.start are "inside":
        // the byte at start belongs to the entry (kept-contiguous rule).
        let mut m = CoordMap::new();
        m.push(entry(0, 1, 1, 100, 50));
        m.apply_indels(&[(100, 1)]);
        assert_eq!(m.entries[0].realised_start_bp, 100);
        assert_eq!(m.entries[0].realised_len_bp, 51);

        let mut m = CoordMap::new();
        m.push(entry(0, 1, 1, 100, 50));
        m.apply_indels(&[(100, -1)]);
        assert_eq!(m.entries[0].realised_start_bp, 100);
        assert_eq!(m.entries[0].realised_len_bp, 49);
    }

    #[test]
    fn shift_after_does_not_extend_right_entry() {
        // F2 regression: an uncovered insertion at pos == entry.start
        // must SHIFT that entry right, not extend it. (Contrast with
        // apply_indels which would extend.)
        let mut m = CoordMap::new();
        m.push(entry(0, 1, 1, 0, 100));
        m.push(entry(0, 1, 2, 100, 100));
        m.shift_after(100, 50);
        assert_eq!(m.entries[0].realised_start_bp, 0);
        assert_eq!(m.entries[0].realised_len_bp, 100);
        assert_eq!(m.entries[1].realised_start_bp, 150);
        assert_eq!(m.entries[1].realised_len_bp, 100); // length unchanged
    }

    #[test]
    fn shift_after_leaves_entries_before_pos_untouched() {
        let mut m = CoordMap::new();
        m.push(entry(0, 1, 1, 0, 100));
        m.push(entry(0, 1, 2, 100, 100));
        m.push(entry(0, 1, 3, 200, 100));
        m.shift_after(200, 50);
        assert_eq!(m.entries[0].realised_start_bp, 0);
        assert_eq!(m.entries[1].realised_start_bp, 100);
        assert_eq!(m.entries[2].realised_start_bp, 250);
    }

    #[test]
    fn apply_indels_tile_stays_contiguous() {
        // Two adjacent entries fully tiling [0, 200). An insertion at the
        // boundary (pos=100) goes to the RIGHT entry; the left entry is
        // unchanged. Sum of lengths must equal the post-insertion span.
        let mut m = CoordMap::new();
        m.push(entry(0, 1, 1, 0, 100));
        m.push(entry(0, 1, 2, 100, 100));
        m.apply_indels(&[(100, 1)]);
        assert_eq!(m.entries[0].realised_start_bp, 0);
        assert_eq!(m.entries[0].realised_len_bp, 100);
        assert_eq!(m.entries[1].realised_start_bp, 100);
        assert_eq!(m.entries[1].realised_len_bp, 101);
        let total: usize = m.entries.iter().map(|e| e.realised_len_bp).sum();
        assert_eq!(total, 201);
    }

    #[test]
    fn apply_indels_at_end_boundary() {
        // insertion at pos == entry.end: unaffected (insertion is AFTER the entry).
        let mut m = CoordMap::new();
        m.push(entry(0, 1, 1, 100, 50));
        m.apply_indels(&[(150, 1)]);
        assert_eq!(m.entries[0].realised_start_bp, 100);
        assert_eq!(m.entries[0].realised_len_bp, 50);

        // deletion at pos == entry.end: unaffected (the byte AT end is not part of entry).
        let mut m = CoordMap::new();
        m.push(entry(0, 1, 1, 100, 50));
        m.apply_indels(&[(150, -1)]);
        assert_eq!(m.entries[0].realised_start_bp, 100);
        assert_eq!(m.entries[0].realised_len_bp, 50);
    }

    #[test]
    fn apply_indels_mixed_multiple_entries() {
        let mut m = CoordMap::new();
        m.push(entry(0, 1, 1, 0, 100));
        m.push(entry(0, 1, 2, 100, 100));
        m.push(entry(0, 2, 1, 200, 100));
        // insertion at pos 50 (inside entry 0); deletion at pos 250 (inside entry 2).
        m.apply_indels(&[(50, 1), (250, -1)]);
        // entry 0: insertion at 50 is inside [0, 100): len 100 → 101
        assert_eq!(m.entries[0].realised_start_bp, 0);
        assert_eq!(m.entries[0].realised_len_bp, 101);
        // entry 1: insertion at 50 < 100 → shift +1; deletion at 250 > 200 → unaffected
        assert_eq!(m.entries[1].realised_start_bp, 101);
        assert_eq!(m.entries[1].realised_len_bp, 100);
        // entry 2: insertion at 50 < 200 → shift +1; deletion at 250 inside [200, 300) → shrink -1
        assert_eq!(m.entries[2].realised_start_bp, 201);
        assert_eq!(m.entries[2].realised_len_bp, 99);
    }
}
