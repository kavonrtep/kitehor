//! Sub-stream seed derivation.
//!
//! Each generation stage gets its own deterministic PRNG seeded as
//! FNV-1a hash of `top_seed.to_le_bytes() || ":" || name`, matching
//! the parent project's convention (see kite2 CLAUDE.md: *"per-case
//! seed is derived deterministically from master_seed:case_id via
//! FNV-1a hash"*).
//!
//! Sub-streams isolate stages: changing the noise rate doesn't
//! perturb template generation, which lets the test suite A/B
//! detector calibrations cleanly.

use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

const FNV_OFFSET: u64 = 0xcbf29ce484222325;
const FNV_PRIME: u64 = 0x100000001b3;

/// FNV-1a hash of `top.to_le_bytes() || ":" || name`.
pub fn derive(top: u64, name: &str) -> u64 {
    let mut h = FNV_OFFSET;
    for b in top.to_le_bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(FNV_PRIME);
    }
    h ^= b':' as u64;
    h = h.wrapping_mul(FNV_PRIME);
    for &b in name.as_bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(FNV_PRIME);
    }
    h
}

pub fn rng(top: u64, name: &str) -> ChaCha20Rng {
    ChaCha20Rng::seed_from_u64(derive(top, name))
}

/// Named sub-streams. Build once with `top_seed`, hand out fresh
/// PRNGs by stage.
#[derive(Debug, Clone, Copy)]
pub struct Streams {
    pub top: u64,
}

impl Streams {
    pub fn new(top: u64) -> Self {
        Self { top }
    }
    pub fn templates(&self) -> ChaCha20Rng {
        rng(self.top, "templates")
    }
    pub fn structure(&self) -> ChaCha20Rng {
        rng(self.top, "structure")
    }
    pub fn wobble(&self) -> ChaCha20Rng {
        rng(self.top, "wobble")
    }
    pub fn events(&self) -> ChaCha20Rng {
        rng(self.top, "events")
    }
    pub fn noise(&self) -> ChaCha20Rng {
        rng(self.top, "noise")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::RngCore;

    #[test]
    fn derive_is_deterministic() {
        assert_eq!(derive(42, "templates"), derive(42, "templates"));
    }

    #[test]
    fn different_names_diverge() {
        assert_ne!(derive(42, "templates"), derive(42, "structure"));
        assert_ne!(derive(42, "templates"), derive(42, "noise"));
    }

    #[test]
    fn different_top_seeds_diverge() {
        assert_ne!(derive(0, "templates"), derive(1, "templates"));
        assert_ne!(derive(1, "templates"), derive(2, "templates"));
    }

    #[test]
    fn streams_isolate_stages() {
        let s = Streams::new(123);
        let mut a = s.templates();
        let mut b = s.noise();
        // First u64 from each stream must differ — they're seeded with
        // different values, so collision is astronomically unlikely.
        assert_ne!(a.next_u64(), b.next_u64());
    }

    #[test]
    fn streams_replay_each_stage() {
        let s = Streams::new(7);
        let v1 = s.templates().next_u64();
        let v2 = s.templates().next_u64();
        assert_eq!(v1, v2, "calling templates() twice must replay the same sequence");
    }
}
