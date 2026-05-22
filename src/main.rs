//! `kitehor` CLI.

use anyhow::{Context, Result};
use clap::Parser;
use kitehor::cli::{
    AnalyzeArgs, Cli, Command, DetectArgs, DetectBatchArgs, HorValidateArgs,
    KitePeriodicityArgs, RuleClassifyArgs, SimulateArgs, SimulateGridArgs, SsrScanArgs,
    SubrepeatScanArgs, SummaryMergeArgs, SynthArgs, SynthBatchArgs, SynthValidateArgs,
};
use kitehor::io::{load_fasta, LoadQc, LoadStatus};
use kitehor::kite::{analyze as kite_analyze, KiteConfig};
use kitehor::simulate::{simulate, SimulateParams};
use log::info;
use rayon::prelude::*;

fn main() -> Result<()> {
    let cli = Cli::parse();
    init_logger(cli.verbose);
    match cli.command {
        Command::KitePeriodicity(args) => run_kite_periodicity(*args),
        Command::Simulate(args) => run_simulate(args),
        Command::SimulateGrid(args) => run_simulate_grid(args),
        Command::SynthValidate(args) => run_synth_validate(args),
        Command::SynthSchema => run_synth_schema(),
        Command::Synth(args) => run_synth(args),
        Command::SynthBatch(args) => run_synth_batch(args),
        Command::Detect(args) => run_detect(args),
        Command::DetectBatch(args) => run_detect_batch(args),
        Command::RuleClassify(args) => run_rule_classify(args),
        Command::SummaryMerge(args) => run_summary_merge(args),
        Command::SsrScan(args) => run_ssr_scan(args),
        Command::SubrepeatScan(args) => run_subrepeat_scan(args),
        Command::HorValidate(args) => run_hor_validate(args),
        Command::Analyze(args) => run_analyze(args),
    }
}

fn run_analyze(args: AnalyzeArgs) -> Result<()> {
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .ok();
    }
    let mut cfg = kitehor::analyze::Config::default();
    cfg.rule.tol = args.rule_tol;
    cfg.rule.min_period = args.rule_min_period;
    cfg.rule.k_max = args.rule_k_max;
    cfg.summary.pure_ssr_pct_threshold = args.pure_ssr_pct_threshold;
    cfg.ssr.ssr_flag_threshold_pct = args.ssr_flag_threshold_pct;
    let report = kitehor::analyze::run(&args.fasta, &args.out, &cfg)?;
    info!(
        "analyze: {} record(s) — hor={} simple_tr={} unresolved={}",
        report.n_records, report.n_hor, report.n_tr, report.n_unresolved
    );
    Ok(())
}

fn run_subrepeat_scan(args: SubrepeatScanArgs) -> Result<()> {
    let cfg = kitehor::subrepeat::Config {
        tol: args.tol,
        window_mult_sub: args.window_mult_sub,
        step_frac: args.step_frac,
        top_n_sub: args.top_n_sub,
        top_n_host: args.top_n_host,
        sub_floor: args.sub_floor,
        window_score_floor: args.window_score_floor,
        min_run: args.min_run,
        host_sub_ratio_min: args.host_sub_ratio_min,
        min_window_bp: args.min_window_bp,
    };
    let n = kitehor::subrepeat::run_subcommand(
        &args.fasta,
        &args.out,
        &args.kite_peaks,
        &cfg,
    )?;
    info!("subrepeat-scan: scanned {n} record(s)");
    Ok(())
}

fn run_hor_validate(args: HorValidateArgs) -> Result<()> {
    let cfg = kitehor::hor_validate::Config {
        period_match_tol: args.period_match_tol,
        min_k_for_density: args.min_k_for_density,
        density_window_tile_frac: args.density_window_tile_frac,
        min_founder_mult: args.min_founder_mult,
        min_density_window_bp: args.min_density_window_bp,
        max_density_windows: args.max_density_windows,
        density_rel_floor: args.density_rel_floor,
        phase_fold_bins: args.phase_fold_bins,
        density_dup_max: args.density_dup_max,
        density_hor_min: args.density_hor_min,
        phase_contrast_dup_min: args.phase_contrast_dup_min,
        phase_contrast_hor_max: args.phase_contrast_hor_max,
        max_tile_bp: args.max_tile_bp,
        min_window_bp: args.min_window_bp,
    };
    let n = kitehor::hor_validate::run_subcommand(
        &args.fasta,
        &args.verdicts,
        &args.global_peaks,
        &args.out,
        &cfg,
    )?;
    info!("hor-validate: wrote {n} row(s)");
    Ok(())
}

fn run_ssr_scan(args: SsrScanArgs) -> Result<()> {
    let specs = kitehor::ssr::parse_motif_min_reps(&args.motif_min_reps)?;
    let cfg = kitehor::ssr::Config {
        ssr_flag_threshold_pct: args.ssr_flag_threshold_pct,
        specs,
        consensus_dimer_copies: args.consensus_dimer_copies,
        consensus_dimer_min_bp: args.consensus_dimer_min_bp,
        consensus_max_monomers: args.consensus_max_monomers,
        consensus_freq_ratio_min: args.consensus_freq_ratio_min,
    };
    let n = kitehor::ssr::run_subcommand(
        &args.fasta,
        &args.out,
        args.kite_peaks.as_deref(),
        &cfg,
    )?;
    info!("ssr-scan: scanned {n} record(s)");
    Ok(())
}

fn run_summary_merge(args: SummaryMergeArgs) -> Result<()> {
    let cfg = kitehor::summary::Config {
        pure_ssr_pct_threshold: args.pure_ssr_pct_threshold,
    };
    let n = kitehor::summary::run_subcommand(
        &args.verdicts,
        &args.subrepeat,
        &args.ssr,
        args.within_tile.as_deref(),
        &args.out,
        &cfg,
    )?;
    info!("summary-merge: wrote {n} row(s)");
    Ok(())
}

fn run_rule_classify(args: RuleClassifyArgs) -> Result<()> {
    let cfg = kitehor::rule_classify::Config {
        tol: args.tol,
        min_period: args.min_period,
        min_cluster_frac: args.min_cluster_frac,
        k_max: args.k_max,
        non_mono_ratio: args.non_mono_ratio,
        founder_floor: args.founder_floor,
        high_k_tile_floor: args.high_k_tile_floor,
        lone_significant_frac: args.lone_significant_frac,
    };
    let n = kitehor::rule_classify::run_subcommand(
        &args.peaks,
        &args.out,
        &cfg,
        args.dump_clusters.as_deref(),
    )?;
    info!("rule-classify: wrote {n} verdict(s)");
    Ok(())
}

// ---------------------------------------------------------------------------
// detect / detect-batch
// ---------------------------------------------------------------------------

fn load_detector_config(path: Option<&std::path::PathBuf>) -> Result<kitehor::detect::DetectorConfig> {
    match path {
        Some(p) => kitehor::detect::DetectorConfig::load(p),
        None => {
            let c = kitehor::detect::DetectorConfig::default();
            c.validate()?;
            Ok(c)
        }
    }
}

fn run_detect(args: DetectArgs) -> Result<()> {
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .ok();
    }
    let cfg = load_detector_config(args.config.as_ref())?;
    // DH9: requesting an export without --viz-dir is a usage error.
    let any_export = args.export_raster || args.export_shift || args.export_edges || args.export_ic;
    if any_export && args.viz_dir.is_none() {
        anyhow::bail!(
            "--export-* flag was supplied without --viz-dir; specify --viz-dir <DIR> for the export root"
        );
    }
    let viz_flags = kitehor::detect::VizFlags {
        viz_dir: args.viz_dir.clone(),
        export_raster: args.export_raster,
        export_shift: args.export_shift,
        export_edges: args.export_edges,
        export_ic: args.export_ic,
    };
    let report = match args.periods.as_ref() {
        Some(p) => kitehor::detect::run_one(
            &args.fasta,
            p,
            &args.out,
            &cfg,
            &viz_flags,
            args.allow_missing_periods,
            args.allow_extra_periods,
        )?,
        None => {
            info!(
                "detect: --periods not supplied; deriving via kite-periodicity \
                 with defaults (writes {:?}.periods.tsv)",
                args.out
            );
            kitehor::detect::run_one_auto(&args.fasta, &args.out, &cfg, &viz_flags)?
        }
    };
    info!(
        "detect: {} array(s), {} segment(s), {} width row(s); prefix {:?}",
        report.n_arrays, report.n_segments, report.n_width_rows, args.out
    );
    Ok(())
}

fn run_detect_batch(args: DetectBatchArgs) -> Result<()> {
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .ok();
    }
    let cfg = load_detector_config(args.config.as_ref())?;
    let any_export = args.export_raster || args.export_shift || args.export_edges || args.export_ic;
    if any_export && args.viz_dir.is_none() {
        anyhow::bail!(
            "--export-* flag was supplied without --viz-dir; specify --viz-dir <DIR> for the export root"
        );
    }
    let viz_flags = kitehor::detect::VizFlags {
        viz_dir: args.viz_dir.clone(),
        export_raster: args.export_raster,
        export_shift: args.export_shift,
        export_edges: args.export_edges,
        export_ic: args.export_ic,
    };
    let n = match args.periods_dir.as_ref() {
        Some(p) => kitehor::detect::run_batch(
            &args.fasta_dir,
            p,
            &args.out_dir,
            &cfg,
            &viz_flags,
            args.allow_missing_periods,
            args.allow_extra_periods,
        )?,
        None => {
            info!(
                "detect-batch: --periods-dir not supplied; deriving periods \
                 per FASTA via kite-periodicity with defaults (writes \
                 <stem>.periods.tsv alongside each output bundle)"
            );
            kitehor::detect::run_batch_auto(&args.fasta_dir, &args.out_dir, &cfg, &viz_flags)?
        }
    };
    info!("detect-batch: processed {n} array(s) into {:?}", args.out_dir);
    Ok(())
}

fn run_synth(args: SynthArgs) -> Result<()> {
    kitehor::synth::run_one(&args.config, &args.out, args.seed, args.diagnostics)?;
    info!(
        "synth: wrote {:?}.fa / .truth.tsv / .periods.tsv",
        args.out
    );
    Ok(())
}

fn run_synth_batch(args: SynthBatchArgs) -> Result<()> {
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .ok();
    }
    let n = kitehor::synth::run_batch(
        &args.config_dir,
        &args.out_dir,
        args.seed_offset,
        args.diagnostics,
    )?;
    info!(
        "synth-batch: wrote {} bundle(s) to {:?}",
        n, args.out_dir
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// synth-validate / synth-schema (M1)
// ---------------------------------------------------------------------------

fn run_synth_validate(args: SynthValidateArgs) -> Result<()> {
    use kitehor::synth::load_and_validate;
    match load_and_validate(&args.config) {
        Ok(cfg) => {
            info!(
                "{:?}: OK — {} block(s), {} template(s), {} event(s)",
                args.config,
                cfg.structure.len(),
                cfg.templates.len(),
                cfg.post_generation.len()
            );
            println!("ok");
            Ok(())
        }
        Err(e) => Err(anyhow::anyhow!("{}", e)),
    }
}

fn run_synth_schema() -> Result<()> {
    use std::io::Write;
    let mut stdout = std::io::stdout().lock();
    stdout.write_all(kitehor::synth::CANONICAL_SCHEMA.as_bytes())?;
    Ok(())
}

fn init_logger(verbosity: u8) {
    let level = match verbosity {
        0 => "warn",
        1 => "info",
        2 => "debug",
        _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level))
        .format_timestamp_secs()
        .init();
}

// ---------------------------------------------------------------------------
// kite-periodicity
// ---------------------------------------------------------------------------

fn run_kite_periodicity(args: KitePeriodicityArgs) -> Result<()> {
    use std::io::Write;

    if args.fasta.is_empty() {
        anyhow::bail!("at least one input FASTA is required");
    }
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .ok();
    }

    let qc = LoadQc {
        min_array_bp: args.qc.min_array_bp,
        max_n_fraction: args.qc.max_n_fraction,
    };
    let cfg = KiteConfig {
        k: args.kmer_size,
        n_bg_replicates: args.bg_replicates,
        score2_threshold: args.score2_threshold,
        min_peak_distance: args.min_peak_distance,
        bg_smoothing_sigma: args.bg_sigma,
    };

    let mut loaded = Vec::new();
    for path in &args.fasta {
        let recs = load_fasta(path, qc).with_context(|| format!("loading {:?}", path))?;
        info!("loaded {} records from {:?}", recs.len(), path);
        loaded.extend(recs);
    }
    let dump_dir = args.dump_profile.clone();
    if let Some(ref d) = dump_dir {
        std::fs::create_dir_all(d)?;
    }

    let ok_records: Vec<&kitehor::sequence::ArrayRecord> = loaded
        .iter()
        .filter_map(|lr| match &lr.status {
            LoadStatus::Ok => Some(&lr.record),
            _ => None,
        })
        .collect();

    let results: Vec<kitehor::kite::KiteResult> = ok_records
        .par_iter()
        .map(|rec| kite_analyze(rec, &cfg, dump_dir.is_some()))
        .collect();

    if let Some(dir) = &dump_dir {
        for r in &results {
            if let (Some(profile), Some(bg)) = (&r.profile, &r.background) {
                let mut p = dir.clone();
                p.push(format!("{}.kite.tsv", r.array_id));
                let mut fh =
                    std::fs::File::create(&p).with_context(|| format!("creating {:?}", &p))?;
                writeln!(fh, "d\tH\tbg")?;
                let n = profile.len().min(bg.len());
                for d in 0..n {
                    if profile[d] > 0.0 || bg[d] > 0.0 {
                        writeln!(fh, "{}\t{:.4}\t{:.4}", d, profile[d], bg[d])?;
                    }
                }
            }
        }
    }

    // --- HOR classification (rule-based, port of rule_proto.py) ---
    let classify_enabled = args.classify;
    let rule_cfg = kitehor::rule_classify::Config::default();
    let rule_verdicts: Vec<kitehor::rule_classify::LegacyVerdict> = if classify_enabled {
        results
            .iter()
            .map(|kr| {
                kitehor::rule_classify::LegacyVerdict::from_verdict(
                    &kitehor::rule_classify::classify(kr, &rule_cfg),
                )
            })
            .collect()
    } else {
        Vec::new()
    };
    if classify_enabled {
        info!(
            "rule classifier: applied to {} record(s)",
            rule_verdicts.len()
        );
    }

    // --- Optional v2-detector periods.tsv emission ---
    if let Some(periods_path) = args.emit_periods.as_ref() {
        let batches: Vec<Vec<kitehor::emit_periods::PeriodsRow>> = if classify_enabled {
            results
                .iter()
                .zip(rule_verdicts.iter())
                .map(|(kr, v)| kitehor::emit_periods::build_rows(kr, Some(v)))
                .collect()
        } else {
            results
                .iter()
                .map(|kr| kitehor::emit_periods::build_rows(kr, None))
                .collect()
        };
        let n = kitehor::emit_periods::write_tsv(periods_path, &batches)?;
        let n_arrays_with_rows = batches.iter().filter(|b| !b.is_empty()).count();
        info!(
            "emit-periods: wrote {} row(s) for {} of {} array(s) to {:?}",
            n,
            n_arrays_with_rows,
            results.len(),
            periods_path
        );
    }

    // Primary TSV.
    if let Some(parent) = args.out.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut out =
        std::fs::File::create(&args.out).with_context(|| format!("creating {:?}", &args.out))?;
    let mut header = String::from(
        "case_id\tarray_length\tn_peaks_kept\
         \tmonomer_size\tscore\
         \tmonomer_size_2\tscore_2\
         \tmonomer_size_3\tscore_3",
    );
    if classify_enabled {
        header.push_str("\tverdict\tfounder\tmultiplicity\ttile\tshare");
    }
    writeln!(out, "{}", header)?;
    let na = "NA";
    for (idx, r) in results.iter().enumerate() {
        let p1 = r.peaks.first();
        let p2 = r.peaks.get(1);
        let p3 = r.peaks.get(2);
        let fmt_period = |p: Option<&kitehor::kite::KitePeak>| {
            p.map(|x| x.period.to_string()).unwrap_or_else(|| na.into())
        };
        let fmt_score = |p: Option<&kitehor::kite::KitePeak>| {
            p.map(|x| format!("{:.10}", x.score))
                .unwrap_or_else(|| na.into())
        };
        let base = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.array_id,
            r.length_bp,
            r.peaks.len(),
            fmt_period(p1),
            fmt_score(p1),
            fmt_period(p2),
            fmt_score(p2),
            fmt_period(p3),
            fmt_score(p3),
        );
        let mut line = base;
        if classify_enabled {
            let rv = rule_verdicts[idx];
            let fmt_opt = |o: Option<usize>| o.map(|x| x.to_string()).unwrap_or_else(|| na.into());
            let fmt_share =
                |s: Option<f64>| s.map(|x| format!("{:.4}", x)).unwrap_or_else(|| na.into());
            line.push_str(&format!(
                "\t{}\t{}\t{}\t{}\t{}",
                rv.as_str(),
                fmt_opt(rv.founder()),
                fmt_opt(rv.multiplicity()),
                fmt_opt(rv.tile()),
                fmt_share(rv.share()),
            ));
        }
        writeln!(out, "{}", line)?;
    }

    // Long-format peaks TSV.
    let out_peaks_path = args.out_peaks.clone().unwrap_or_else(|| {
        let mut p = args.out.as_os_str().to_owned();
        p.push(".peaks.tsv");
        std::path::PathBuf::from(p)
    });
    let mut out_peaks = std::fs::File::create(&out_peaks_path)
        .with_context(|| format!("creating {:?}", &out_peaks_path))?;
    writeln!(
        out_peaks,
        "case_id\tarray_length\trank\tperiod\tpeak_height\tscore\tscore2\tscore2_norm\tbackground"
    )?;
    for r in &results {
        for (i, p) in r.peaks.iter().enumerate() {
            writeln!(
                out_peaks,
                "{}\t{}\t{}\t{}\t{:.4}\t{:.10}\t{:.10}\t{:.10}\t{:.4}",
                r.array_id,
                r.length_bp,
                i + 1,
                p.period,
                p.peak_height,
                p.score,
                p.score2,
                p.score2_norm,
                p.background
            )?;
        }
    }

    info!(
        "wrote {} record(s) to {:?}; long-format to {:?}",
        results.len(),
        args.out,
        out_peaks_path
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// simulate
// ---------------------------------------------------------------------------

fn run_simulate(args: SimulateArgs) -> Result<()> {
    use std::io::Write;
    let params = SimulateParams {
        monomer_len: args.monomer_size,
        hor_order: args.multiplicity,
        n_blocks: args.copies,
        sub_rate_intra: args.sub_rate_intra,
        sub_rate_inter: args.sub_rate_inter,
        submono_k: args.submono_k,
        seed: args.seed,
        ..SimulateParams::default()
    };
    let (arr, truth, _monomers, _events) =
        simulate(&args.case_id, &params).map_err(|e| anyhow::anyhow!(e.to_string()))?;
    if let Some(parent) = args.out.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    {
        let mut f = std::fs::File::create(&args.out)?;
        writeln!(f, ">{}", arr.id)?;
        for chunk in arr.seq.chunks(60) {
            f.write_all(chunk)?;
            f.write_all(b"\n")?;
        }
    }
    let truth_path = {
        let mut p = args.out.as_os_str().to_owned();
        p.push(".truth.tsv");
        std::path::PathBuf::from(p)
    };
    {
        let mut f = std::fs::File::create(&truth_path)?;
        writeln!(
            f,
            "case_id\tmonomer_len\thor_order\tn_blocks\tsub_rate_intra\tsub_rate_inter\tsubmono_k\tseed\tarray_length\tn_monomers"
        )?;
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            truth.case_id,
            truth.monomer_len,
            truth.hor_order,
            truth.n_blocks,
            truth.sub_rate_intra,
            truth.sub_rate_inter,
            truth.submono_k,
            truth.seed,
            truth.array_length,
            truth.n_monomers,
        )?;
    }
    info!(
        "wrote {} bp to {:?} and truth to {:?}",
        arr.length, args.out, truth_path
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// simulate-grid
// ---------------------------------------------------------------------------

fn run_simulate_grid(args: SimulateGridArgs) -> Result<()> {
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .ok();
    }
    kitehor::simulate_grid::run_grid(&args.params, &args.outdir, args.seed)
}
