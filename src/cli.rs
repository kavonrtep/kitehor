//! CLI argument parsing (clap derive).

use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(
    name = "kitehor",
    version,
    about = "Kite-first probabilistic HOR detector"
)]
pub struct Cli {
    /// Verbosity level: pass once for INFO, twice for DEBUG, thrice for TRACE.
    #[arg(short, long, action = clap::ArgAction::Count, global = true)]
    pub verbose: u8,

    #[command(subcommand)]
    pub command: Command,
}

#[derive(Debug, Subcommand)]
pub enum Command {
    /// Kite periodicity scan + (optional) HOR classifier.
    /// Primary entry point of the tool.
    KitePeriodicity(Box<KitePeriodicityArgs>),
    /// Generate one synthetic HOR / tandem-repeat array (for testing
    /// and training-set construction).
    Simulate(SimulateArgs),
    /// Generate a grid of synthetic arrays from a params TSV.
    /// Writes sequences.fasta, truth.tsv, monomers.tsv, events.tsv,
    /// and alternatives.tsv to the output directory.
    SimulateGrid(SimulateGridArgs),
    /// Validate a v2 YAML simulator config against the canonical
    /// schema (`docs/new/simulator_schema.json`) and MVP invariants.
    /// Exits non-zero on first validation error.
    SynthValidate(SynthValidateArgs),
    /// Print the canonical JSON Schema to stdout.
    SynthSchema,
    /// Generate one synthetic tandem-repeat array from a YAML config.
    /// Writes {PREFIX}.fa, {PREFIX}.truth.tsv, {PREFIX}.periods.tsv.
    Synth(SynthArgs),
    /// Run `synth` over every `*.yaml` in a directory (parallel).
    /// `.deferred.yaml` placeholders are skipped.
    SynthBatch(SynthBatchArgs),
    /// v2 line-width detector. Reads FASTA + period candidates, emits
    /// {PREFIX}.properties.tsv, .segments.tsv, .width_features.tsv,
    /// .consensus.fa, and (optionally) per-array visualisation TSVs/PNGs.
    /// Experimental: per-category benchmark accuracy is still being
    /// calibrated — see docs/new/detect_impl_plan.md.
    Detect(DetectArgs),
    /// Run `detect` over every `<stem>.fa` in a directory paired with
    /// `<stem>.periods.tsv` (parallel).
    DetectBatch(DetectBatchArgs),
    /// Rule-based HOR/simple_tr/unresolved classifier (port of
    /// `tools/rule_proto/rule_proto.py`). Reads a kite peaks
    /// long-format TSV and emits `<prefix>.verdicts.tsv` with one row
    /// per record.
    RuleClassify(RuleClassifyArgs),
    /// Outer-join rule-classify verdicts + subrepeat-scan + ssr-scan
    /// (optionally + hor-validate) into a single summary TSV with
    /// `combined_class`. Port of `tools/rule_proto/summary.py`.
    SummaryMerge(SummaryMergeArgs),
    /// Short-motif tandem repeat (SSR) scanner. Port of
    /// `tools/rule_proto/ssr_scan.py`. Emits `<prefix>.ssr.tsv` and
    /// `<prefix>.ssr.regions.tsv`.
    SsrScan(SsrScanArgs),
    /// Unified spatial-localization subrepeat detector (port of
    /// `tools/rule_proto/tandem_validate.py`, spec v5). Replaces the
    /// older `subrepeat-scan` + `hor-validate` stages in the rule-proto
    /// pipeline. Emits `<prefix>.tandem_validate.tsv`.
    TandemValidate(TandemValidateArgs),
    /// End-to-end pipeline. Runs kite-periodicity → rule-classify →
    /// (tandem-validate ‖ ssr-scan) → summary-merge on one FASTA.
    /// Always emits every per-stage TSV under `<prefix>.*.tsv`.
    Analyze(AnalyzeArgs),
}

// ---------------------------------------------------------------------------
// Shared QC flags
// ---------------------------------------------------------------------------

#[derive(Debug, Args, Clone)]
pub struct QcOpts {
    /// Minimum array length (bp); shorter records are flagged.
    #[arg(long, default_value_t = 5_000)]
    pub min_array_bp: usize,

    /// Maximum fraction of Ns; records above are flagged.
    #[arg(long, default_value_t = 0.20)]
    pub max_n_fraction: f64,
}

impl Default for QcOpts {
    fn default() -> Self {
        Self {
            min_array_bp: 5_000,
            max_n_fraction: 0.20,
        }
    }
}

// ---------------------------------------------------------------------------
// kite-periodicity
// ---------------------------------------------------------------------------

#[derive(Debug, Args)]
pub struct KitePeriodicityArgs {
    /// Input FASTA file(s). For TideCluster-style dimers, pre-trim to
    /// the first half externally.
    pub fasta: Vec<PathBuf>,

    /// Primary output TSV.
    #[arg(short, long)]
    pub out: PathBuf,

    /// k-mer size for the kite histogram. kite.R default = 6.
    #[arg(short = 'k', long, default_value_t = 6)]
    pub kmer_size: usize,

    /// Number of composition-matched random sequences for the noise
    /// envelope. kite.R default = 10.
    #[arg(short = 'N', long = "bg-replicates", default_value_t = 10)]
    pub bg_replicates: usize,

    /// Score2 threshold (kite.R's `threshold` argument).
    #[arg(long, default_value_t = 0.001)]
    pub score2_threshold: f64,

    /// Minimum peak separation (bp).
    #[arg(long, default_value_t = 1)]
    pub min_peak_distance: usize,

    /// Override the gaussian smoothing sigma used to approximate
    /// `smooth.spline` on the random-background envelope. Default = 10.
    #[arg(long)]
    pub bg_sigma: Option<f64>,

    /// Optional: write a FASTA-like periodogram bundle to this path.
    /// Each input record produces two lines: a `>case_id|H` header
    /// followed by the dense neighbour-distance histogram (one integer
    /// per period 1..N), and a `>case_id|bg` header followed by the
    /// composition-matched random background envelope (one float per
    /// period). Mirrors the data shape TideCluster keeps in its
    /// in-memory `profile_list` (see `tarean/kite.R`); useful for
    /// plotting periodograms with any text-aware tool. Off by default.
    #[arg(long, value_name = "PATH")]
    pub periodogram: Option<PathBuf>,

    /// Optional: long-format peaks TSV (one row per kept peak per
    /// record). Defaults to `<out>.peaks.tsv` when omitted.
    #[arg(long)]
    pub out_peaks: Option<PathBuf>,

    /// Optional: emit a v2 `periods.tsv` consumable by
    /// `kitehor detect --periods`. Score mapping (see
    /// `src/emit_periods.rs`): with `--classify`, founder → 0.95,
    /// tile → 0.90, other top-3 peaks → 0.60; ambiguous verdicts
    /// (`Unresolved`) emit top-3 at 0.50 / 0.40 / 0.30. Without
    /// `--classify`, raw kite top-3 peaks emit at 0.60. NoSignal
    /// produces no rows.
    #[arg(long, value_name = "PATH")]
    pub emit_periods: Option<PathBuf>,

    /// Apply the rule-based HOR classifier (port of
    /// `tools/rule_proto/rule_proto.py`). Adds per-record columns
    /// `verdict`, `founder`, `multiplicity`, `tile`, `share`.
    #[arg(long)]
    pub classify: bool,

    #[command(flatten)]
    pub qc: QcOpts,

    /// Number of rayon worker threads (0 = auto).
    #[arg(long, default_value_t = 0)]
    pub threads: usize,
}

// ---------------------------------------------------------------------------
// simulate / simulate-grid
// ---------------------------------------------------------------------------

#[derive(Debug, Args)]
pub struct SimulateArgs {
    /// Base monomer length (bp).
    #[arg(long, default_value_t = 171)]
    pub monomer_size: usize,

    /// HOR multiplicity (1 = plain tandem repeat).
    #[arg(long, default_value_t = 12)]
    pub multiplicity: usize,

    /// Number of HOR copies in the array.
    #[arg(long, default_value_t = 100)]
    pub copies: usize,

    /// Per-base substitution rate within each monomer copy.
    #[arg(long, default_value_t = 0.05)]
    pub sub_rate_intra: f64,

    /// Per-base substitution rate between HOR-position founders.
    #[arg(long, default_value_t = 0.03)]
    pub sub_rate_inter: f64,

    /// Sub-motif tiling factor (1 = no internal sub-period).
    #[arg(long, default_value_t = 1)]
    pub submono_k: usize,

    /// Seed for the simulator's PRNG. Any change → different sequence.
    #[arg(long, default_value_t = 0)]
    pub seed: u64,

    /// Record id placed in the FASTA header.
    #[arg(long, default_value = "sim_0000")]
    pub case_id: String,

    /// Output FASTA path. Truth metadata is written next to it with
    /// `.truth.tsv` appended.
    #[arg(short, long)]
    pub out: PathBuf,
}

#[derive(Debug, Args)]
pub struct SimulateGridArgs {
    /// Parameter TSV. Columns: case_id, monomer_len, hor_order, n_blocks,
    /// sub_rate_intra, sub_rate_inter, indel_rate_intra, indel_rate_inter,
    /// block_conversions, monomer_conversions, submono_k, seed.
    /// Same schema as `ground_truth/params.tsv`.
    #[arg(short, long)]
    pub params: PathBuf,

    /// Output directory. Writes sequences.fasta, truth.tsv,
    /// monomers.tsv, events.tsv, alternatives.tsv here.
    #[arg(short, long)]
    pub outdir: PathBuf,

    /// Master seed used when a row's `seed` column is blank. The
    /// per-case seed is then derived as FNV-1a hash of "master:case_id".
    #[arg(long, default_value_t = 42)]
    pub seed: u64,

    /// Number of rayon worker threads (0 = let rayon decide).
    #[arg(long, default_value_t = 0)]
    pub threads: usize,
}

// ---------------------------------------------------------------------------
// synth-* (v2 simulator)
// ---------------------------------------------------------------------------

#[derive(Debug, Args)]
pub struct SynthValidateArgs {
    /// YAML config to validate.
    pub config: PathBuf,
}

#[derive(Debug, Args)]
pub struct SynthArgs {
    /// YAML config file.
    pub config: PathBuf,
    /// Output prefix (REQUIRED). Writes PREFIX.fa, PREFIX.truth.tsv,
    /// PREFIX.periods.tsv. There is no fallback to `global.output` in
    /// the YAML (silently ignored — see plan §0 A3).
    #[arg(short, long, required = true)]
    pub out: PathBuf,
    /// Override the YAML's `seed`.
    #[arg(long)]
    pub seed: Option<u64>,
    /// Also emit PREFIX.diagnostics.json.
    #[arg(long)]
    pub diagnostics: bool,
}

// ---------------------------------------------------------------------------
// detect / detect-batch (v2 line-width detector)
// ---------------------------------------------------------------------------

#[derive(Debug, Args)]
pub struct DetectArgs {
    /// Input FASTA (one or many records).
    pub fasta: PathBuf,
    /// Period candidates TSV (matches `kitehor synth` output schema).
    /// When omitted, kite-periodicity runs internally with default
    /// settings to derive period candidates; the derived periods are
    /// persisted to `<out>.periods.tsv` as part of the output bundle
    /// and the implicit equivalent of `--allow-missing-periods` is
    /// applied (records that fail kite QC end up classified as
    /// `ambiguous`).
    #[arg(long)]
    pub periods: Option<PathBuf>,
    /// Output prefix; writes PREFIX.properties.tsv, .segments.tsv,
    /// .width_features.tsv, .consensus.fa, and .diagnostics.json.
    #[arg(short, long, required = true)]
    pub out: PathBuf,
    /// Override default DetectorConfig via TOML.
    #[arg(long)]
    pub config: Option<PathBuf>,
    /// Per-array matrix/signal export directory.
    #[arg(long)]
    pub viz_dir: Option<PathBuf>,
    /// Granular viz flags. Require `--viz-dir` (DH9).
    #[arg(long)]
    pub export_raster: bool,
    #[arg(long)]
    pub export_shift: bool,
    #[arg(long)]
    pub export_edges: bool,
    #[arg(long)]
    pub export_ic: bool,
    /// Downgrade missing-periods errors to warnings (DH5).
    #[arg(long)]
    pub allow_missing_periods: bool,
    /// Downgrade "periods TSV has rows whose array_id matches no
    /// FASTA record" errors to warnings (Review-2026-05-16 #11).
    #[arg(long)]
    pub allow_extra_periods: bool,
    /// Number of rayon worker threads (0 = auto).
    #[arg(long, default_value_t = 0)]
    pub threads: usize,
}

#[derive(Debug, Args)]
pub struct DetectBatchArgs {
    /// Directory of FASTA files (`<stem>.fa`).
    #[arg(long)]
    pub fasta_dir: PathBuf,
    /// Directory of period TSVs (`<stem>.periods.tsv`). When omitted,
    /// kite-periodicity runs internally per FASTA to derive period
    /// candidates and the derived periods are written to
    /// `<out_dir>/<stem>.periods.tsv` alongside the detector outputs;
    /// `--allow-missing-periods` / `--allow-extra-periods` are then
    /// ignored because the auto path generates exactly the rows it
    /// uses.
    #[arg(long)]
    pub periods_dir: Option<PathBuf>,
    /// Output directory. Each `<stem>` produces `<stem>.{properties,segments,width_features}.tsv`.
    #[arg(long)]
    pub out_dir: PathBuf,
    /// Override default DetectorConfig via TOML.
    #[arg(long)]
    pub config: Option<PathBuf>,
    /// Per-array matrix/signal export directory.
    #[arg(long)]
    pub viz_dir: Option<PathBuf>,
    /// Granular viz flags (DH10 — mirror `detect`).
    #[arg(long)]
    pub export_raster: bool,
    #[arg(long)]
    pub export_shift: bool,
    #[arg(long)]
    pub export_edges: bool,
    #[arg(long)]
    pub export_ic: bool,
    /// Downgrade missing-periods errors to warnings (DH5).
    #[arg(long)]
    pub allow_missing_periods: bool,
    /// Allow periods TSVs without a matching FASTA stem in the batch
    /// (DH11 — symmetric pairing is the default).
    #[arg(long)]
    pub allow_extra_periods: bool,
    #[arg(long, default_value_t = 0)]
    pub threads: usize,
}

// ---------------------------------------------------------------------------
// rule-classify
// ---------------------------------------------------------------------------

#[derive(Debug, Args)]
pub struct RuleClassifyArgs {
    /// Kite peaks long-format TSV (one row per kept peak per record).
    /// Required columns: `case_id, rank, period, score2_norm`.
    pub peaks: PathBuf,

    /// Output prefix. Writes `<prefix>.verdicts.tsv`.
    #[arg(short, long)]
    pub out: PathBuf,

    /// Cluster-gap + multiplicity rounding tolerance (relative).
    #[arg(long, default_value_t = 0.015)]
    pub tol: f64,

    /// Drop peaks with `period < min_period` (k-mer floor).
    #[arg(long, default_value_t = 20)]
    pub min_period: usize,

    /// Drop clusters with `total_score < min_cluster_frac × max_cluster_score`.
    #[arg(long, default_value_t = 0.01)]
    pub min_cluster_frac: f64,

    /// Maximum multiplicity considered.
    #[arg(long, default_value_t = 30)]
    pub k_max: usize,

    /// Case B k=2: a longer `k×top` cluster qualifies as tile when its
    /// `total_score / top.total_score` exceeds this.
    #[arg(long, default_value_t = 0.5)]
    pub non_mono_ratio: f64,

    /// Case A: founder cluster score must be >= `founder_floor × top.score`.
    #[arg(long, default_value_t = 0.1)]
    pub founder_floor: f64,

    /// Case B k>=3: tile cluster score floor (× top.score).
    #[arg(long, default_value_t = 0.05)]
    pub high_k_tile_floor: f64,

    /// Lone-significant-cluster fallback floor (× top.score). After
    /// Case A and Case B fail AND no harmonic multiples exist, if
    /// exactly one cluster passes this floor, call the array a
    /// simple_tr with `reason = lone_significant_cluster`.
    #[arg(long, default_value_t = 0.1)]
    pub lone_significant_frac: f64,

    /// Minimum number of tile copies (`array_length / tile`) required
    /// to lock in a `hor` verdict. Calls with fewer copies are demoted
    /// to `unresolved` with `reason = insufficient_tile_copies:n=<f>:<min>`.
    /// 0 disables the gate.
    #[arg(long, default_value_t = 6)]
    pub min_tile_copies: usize,

    /// Optional: write per-case cluster dumps to this directory.
    #[arg(long, value_name = "DIR")]
    pub dump_clusters: Option<PathBuf>,
}

// ---------------------------------------------------------------------------
// analyze (end-to-end orchestrator)
// ---------------------------------------------------------------------------

#[derive(Debug, Args)]
pub struct AnalyzeArgs {
    /// Input FASTA.
    pub fasta: PathBuf,
    /// Output prefix. Writes 6 TSVs: `<prefix>.kite.tsv`,
    /// `.kite.peaks.tsv`, `.verdicts.tsv`, `.tandem_validate.tsv`,
    /// `.ssr.tsv`, `.ssr.regions.tsv`, `.summary.tsv`.
    #[arg(short, long, required = true)]
    pub out: PathBuf,
    /// Number of rayon worker threads (0 = auto).
    #[arg(long, default_value_t = 0)]
    pub threads: usize,
    // Stage-prefixed tunables.
    #[arg(long, default_value_t = 0.015)]
    pub rule_tol: f64,
    #[arg(long, default_value_t = 20)]
    pub rule_min_period: usize,
    #[arg(long, default_value_t = 30)]
    pub rule_k_max: usize,
    /// Minimum number of tile copies (`array_length / tile`) required
    /// to lock in a `hor` verdict. 0 disables.
    #[arg(long, default_value_t = 6)]
    pub rule_min_tile_copies: usize,
    /// `pure_ssr` fires when `ssr_raw_total_coverage_pct ≥` this.
    #[arg(long, default_value_t = 80.0)]
    pub pure_ssr_pct_threshold: f64,
    /// `<base>_with_ssr` partner classes fire when
    /// `ssr_raw_total_coverage_pct ≥` this (but below the pure
    /// threshold). v0.11+.
    #[arg(long, default_value_t = 30.0)]
    pub ssr_has_pct_threshold: f64,
    /// Per-record `ssr_flag = yes` when raw total SSR coverage of the
    /// array `≥` this. Default matches `ssr_has_pct_threshold` so
    /// the recomputed flag column tracks the cascade.
    #[arg(long, default_value_t = 30.0)]
    pub ssr_flag_threshold_pct: f64,

    /// Optional: write a FASTA-like periodogram bundle alongside the
    /// per-stage TSVs (same format as `kitehor kite-periodicity
    /// --periodogram`). Off by default.
    #[arg(long, value_name = "PATH")]
    pub periodogram: Option<PathBuf>,
}

// ---------------------------------------------------------------------------
// tandem-validate
// ---------------------------------------------------------------------------

#[derive(Debug, Args)]
pub struct TandemValidateArgs {
    /// Input FASTA.
    pub fasta: PathBuf,
    /// rule-classify verdicts TSV.
    #[arg(long, required = true)]
    pub verdicts: PathBuf,
    /// Global kite peaks TSV (long-format).
    #[arg(long, required = true)]
    pub peaks: PathBuf,
    /// Output prefix. Writes `<prefix>.tandem_validate.tsv`.
    #[arg(short, long, required = true)]
    pub out: PathBuf,

    // Candidate selection
    #[arg(long, default_value_t = 20)]
    pub cand_min_period: usize,
    #[arg(long, default_value_t = 0.0)]
    pub cand_score_floor: f64,
    #[arg(long, default_value_t = 0.03)]
    pub cand_rel_score_floor: f64,
    #[arg(long, default_value_t = 5)]
    pub cand_top_n: usize,
    #[arg(long, default_value_t = 1.0 / 3.0)]
    pub host_inside_ratio: f64,
    #[arg(long, default_value_t = 0.05)]
    pub founder_tol: f64,

    // Window sizing
    #[arg(long, default_value_t = 1.0 / 3.0)]
    pub window_host_frac: f64,
    #[arg(long, default_value_t = 3.0)]
    pub window_cand_mult: f64,
    #[arg(long, default_value_t = 200)]
    pub min_window_bp: usize,

    // Presence
    #[arg(long, default_value_t = 0.02)]
    pub period_match_tol: f64,
    #[arg(long, default_value_t = 0.3)]
    pub window_score_floor: f64,
    #[arg(long, default_value_t = 0.2)]
    pub presence_rel_floor: f64,

    // Binning + decision thresholds
    #[arg(long, default_value_t = 10)]
    pub n_bins: usize,
    #[arg(long, default_value_t = 0.35)]
    pub density_dup_max: f64,
    #[arg(long, default_value_t = 0.7)]
    pub density_hor_min: f64,
    #[arg(long, default_value_t = 0.4)]
    pub contrast_dup_min: f64,
    #[arg(long, default_value_t = 0.15)]
    pub contrast_hor_max: f64,
    #[arg(long, default_value_t = 3)]
    pub min_present_windows: usize,

    /// Number of rayon worker threads (0 = auto).
    #[arg(long, default_value_t = 0)]
    pub threads: usize,
}

// ---------------------------------------------------------------------------
// ssr-scan
// ---------------------------------------------------------------------------

#[derive(Debug, Args)]
pub struct SsrScanArgs {
    /// Input FASTA.
    pub fasta: PathBuf,
    /// Output prefix. Writes `<prefix>.ssr.tsv` and `<prefix>.ssr.regions.tsv`.
    #[arg(short, long, required = true)]
    pub out: PathBuf,
    /// Optional kite peaks TSV. When supplied, the authoritative scan
    /// uses the consensus-dimer path; otherwise it falls back to raw.
    #[arg(long)]
    pub kite_peaks: Option<PathBuf>,
    /// `ssr_flag = yes` when `dominant_motif_coverage_pct ≥` this.
    #[arg(long, default_value_t = 30.0)]
    pub ssr_flag_threshold_pct: f64,
    /// `monomer × n` copies for the consensus dimer.
    #[arg(long, default_value_t = 4)]
    pub consensus_dimer_copies: usize,
    /// Minimum consensus dimer length (bp) (extended by repeating).
    #[arg(long, default_value_t = 30)]
    pub consensus_dimer_min_bp: usize,
    /// Top-K canonical-distinct consensus monomers considered per record.
    #[arg(long, default_value_t = 3)]
    pub consensus_max_monomers: usize,
    /// Stop extracting consensus monomers when count drops below this
    /// fraction of the top k-mer's count.
    #[arg(long, default_value_t = 0.3)]
    pub consensus_freq_ratio_min: f64,
    /// Per-motif-length minimum repeat counts as `"L:reps,L:reps,…"`.
    /// Default: TideCluster — `"1:20,2:9,3:6,4:5,5:5,6:5,7:5,8:5,9:5,10:5,11:5,12:5,13:5,14:5"`.
    #[arg(
        long,
        default_value = "1:20,2:9,3:6,4:5,5:5,6:5,7:5,8:5,9:5,10:5,11:5,12:5,13:5,14:5"
    )]
    pub motif_min_reps: String,
}

#[derive(Debug, Args)]
pub struct SummaryMergeArgs {
    /// rule-classify verdicts TSV (`case_id` column).
    #[arg(long, required = true)]
    pub verdicts: PathBuf,
    /// tandem-validate TSV (`record_id` column).
    #[arg(long, required = true)]
    pub tandem_validate: PathBuf,
    /// ssr-scan TSV (`record_id` column).
    #[arg(long, required = true)]
    pub ssr: PathBuf,
    /// Output prefix. Writes `<prefix>.summary.tsv`.
    #[arg(short, long, required = true)]
    pub out: PathBuf,
    /// `pure_ssr` fires when `ssr_raw_total_coverage_pct ≥` this.
    #[arg(long, default_value_t = 80.0)]
    pub pure_ssr_pct_threshold: f64,
    /// `<base>_with_ssr` partner classes (`hor_with_ssr`,
    /// `tr_with_ssr`, `unresolved_with_ssr`, `tr_with_subrepeat_with_ssr`)
    /// fire when `ssr_raw_total_coverage_pct ≥` this (but below the
    /// pure threshold). v0.11+.
    #[arg(long, default_value_t = 30.0)]
    pub ssr_has_pct_threshold: f64,
}

#[derive(Debug, Args)]
pub struct SynthBatchArgs {
    /// Directory of `*.yaml` configs (`.deferred.yaml` skipped).
    #[arg(long)]
    pub config_dir: PathBuf,
    /// Output directory. Each config writes `<stem>.fa /
    /// .truth.tsv / .periods.tsv` here.
    #[arg(long)]
    pub out_dir: PathBuf,
    /// Per-config seed = FNV-1a(seed_offset, filename). A fixed
    /// `seed_offset` keeps the corpus byte-reproducible across runs.
    #[arg(long, default_value_t = 0)]
    pub seed_offset: u64,
    /// Also emit `<stem>.diagnostics.json` per config.
    #[arg(long)]
    pub diagnostics: bool,
    /// Number of rayon worker threads (0 = auto).
    #[arg(long, default_value_t = 0)]
    pub threads: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_help_renders() {
        Cli::command().debug_assert();
    }

    #[test]
    fn classify_parses() {
        let argv = [
            "kitehor",
            "kite-periodicity",
            "x.fa",
            "-o",
            "out.tsv",
            "--classify",
        ];
        let cli = Cli::try_parse_from(argv).unwrap();
        match cli.command {
            Command::KitePeriodicity(args) => {
                assert!(args.classify);
                assert_eq!(args.fasta.len(), 1);
            }
            _ => panic!("expected KitePeriodicity"),
        }
    }
}
