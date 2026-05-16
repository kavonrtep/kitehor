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

    /// Optional: write per-record `H[d]` + `bg[d]` profiles to this
    /// directory as `<case_id>.kite.tsv` (large output; off by default).
    #[arg(long)]
    pub dump_profile: Option<PathBuf>,

    /// Optional: long-format peaks TSV (one row per kept peak per
    /// record). Defaults to `<out>.peaks.tsv` when omitted.
    #[arg(long)]
    pub out_peaks: Option<PathBuf>,

    // -- Rule-based HOR layer (legacy, hor_call.rs) --------------------
    /// Disable the rule-based HOR layer (kite-peaks-only output).
    #[arg(long)]
    pub no_hor_call: bool,

    /// Rule layer: max multiplicity explored.
    #[arg(long, default_value_t = 30)]
    pub hor_qmax: usize,

    /// Rule layer: minimum peaks in a family for the HOR verdict.
    #[arg(long, default_value_t = 3)]
    pub hor_min_family_size: usize,

    /// Rule layer: family_score / total_score floor for the HOR verdict.
    #[arg(long, default_value_t = 0.50)]
    pub hor_min_family_share: f64,

    /// Rule layer: top-peak / second-peak score ratio for the tandem
    /// verdict.
    #[arg(long, default_value_t = 3.0)]
    pub hor_dominance: f64,

    /// Rule layer: relative band (± × top1 period) for the jitter
    /// detector that flags variable-length tandems as `unresolved`.
    #[arg(long, default_value_t = 0.15)]
    pub hor_jitter_tol: f64,

    /// Rule layer: # of top peaks within the jitter band that triggers
    /// `unresolved`.
    #[arg(long, default_value_t = 4)]
    pub hor_jitter_thr: usize,

    /// Rule layer: minimum `tile_score / founder_score` for the HOR
    /// verdict.
    #[arg(long, default_value_t = 0.15)]
    pub hor_min_tile_founder_ratio: f64,

    // -- HOR classification ---------------------------------------
    /// Apply the HOR classifier. By default this is the rule-based
    /// classifier (`src/rule.rs`): d1 = k × p_n for some k ∈ [2,
    /// qmax], with p_n a top-N kite peak. Adds per-record columns
    /// `verdict`, `founder`, `multiplicity`, `tile`, `share`.
    /// Pass `--use-ml-classifier` to fall back to the legacy
    /// probabilistic pipeline.
    #[arg(long)]
    pub classify: bool,

    /// Use the legacy ML classifier (random forest + Platt + verdict
    /// logic from earlier kitehor versions) instead of the rule.
    /// Adds the ML-specific columns: `hor_score`, `hor_score_raw`,
    /// `k_pred`, `recovered`, `h_d1`, `h_founder`. Requires the
    /// random-forest JSON artefacts shipped under `models/`.
    #[arg(long)]
    pub use_ml_classifier: bool,

    /// Rule layer: maximum HOR multiplicity considered. Default 30.
    #[arg(long, default_value_t = 30)]
    pub rule_qmax: usize,

    /// Rule layer: founder candidate must be among the top-N kite
    /// peaks by score. Default 3 (the user-validated value).
    #[arg(long, default_value_t = 3)]
    pub rule_top_n: usize,

    /// Supplementary HOR-coverage QC (rule path only). For each HOR
    /// call, slide a tile-length window across the array with step =
    /// tile and compute Levenshtein identity vs the first tile. Adds
    /// columns: `cov_mean`, `cov_pass_70/80/90`, `cov_first_half`,
    /// `cov_second_half`, `cov_min`, `cov_max`, `cov_n_tiles`. Records
    /// not called HOR get NA. Adds modest runtime (~ms per record for
    /// typical tile sizes; longer for kb-scale tiles).
    #[arg(long)]
    pub coverage: bool,

    /// ML override: classifier config TOML.
    #[arg(long, value_name = "PATH")]
    pub classifier_config: Option<PathBuf>,

    /// ML override: HOR-score RF model JSON.
    #[arg(long, value_name = "PATH")]
    pub hor_model: Option<PathBuf>,

    /// ML override: k-predictor RF model JSON.
    #[arg(long, value_name = "PATH")]
    pub k_model: Option<PathBuf>,

    /// ML option: skip the homology probe (`h_d1`, `h_founder`); the
    /// model falls back to the training-set imputation medians.
    /// Only has effect under `--use-ml-classifier`.
    #[arg(long)]
    pub no_homology: bool,

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
    #[arg(long)]
    pub periods: PathBuf,
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
    /// Number of rayon worker threads (0 = auto).
    #[arg(long, default_value_t = 0)]
    pub threads: usize,
}

#[derive(Debug, Args)]
pub struct DetectBatchArgs {
    /// Directory of FASTA files (`<stem>.fa`).
    #[arg(long)]
    pub fasta_dir: PathBuf,
    /// Directory of period TSVs (`<stem>.periods.tsv`).
    #[arg(long)]
    pub periods_dir: PathBuf,
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
