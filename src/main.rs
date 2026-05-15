//! `kitehor` CLI.

use anyhow::{Context, Result};
use clap::Parser;
use kitehor::cli::{Cli, Command, KitePeriodicityArgs, SimulateArgs, SimulateGridArgs};
use kitehor::classifier::{ClassifierConfig, RandomForest};
use kitehor::classify::{classify as run_classify, Verdict};
use kitehor::features::{build_features, FeatureRow};
use kitehor::hor_call::{classify as hor_classify, HorCallConfig};
use kitehor::io::{load_fasta, LoadQc, LoadStatus};
use kitehor::kite::{analyze as kite_analyze, KiteConfig};
use kitehor::monomer_model::{probe_period, MonomerModelConfig};
use kitehor::simulate::{simulate, SimulateParams};
use log::{info, warn};
use rayon::prelude::*;

fn main() -> Result<()> {
    let cli = Cli::parse();
    init_logger(cli.verbose);
    match cli.command {
        Command::KitePeriodicity(args) => run_kite_periodicity(args),
        Command::Simulate(args) => run_simulate(args),
        Command::SimulateGrid(args) => run_simulate_grid(args),
    }
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
    let hor_cfg = HorCallConfig {
        qmax: args.hor_qmax,
        min_family_size: args.hor_min_family_size,
        min_family_share: args.hor_min_family_share,
        dominance: args.hor_dominance,
        jitter_tol: args.hor_jitter_tol,
        jitter_thr: args.hor_jitter_thr,
        min_tile_founder_ratio: args.hor_min_tile_founder_ratio,
        ..HorCallConfig::default()
    };
    let hor_enabled = !args.no_hor_call;

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
                let mut fh = std::fs::File::create(&p)
                    .with_context(|| format!("creating {:?}", &p))?;
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

    // --- Optional probabilistic classifier (kite-first RF + Platt) ---
    let classify_enabled = args.classify;
    let verdicts: Vec<(FeatureRow, Verdict)> = if classify_enabled {
        let cls_cfg = match &args.classifier_config {
            Some(p) => ClassifierConfig::load(p)?,
            None => ClassifierConfig::default_baked()?,
        };
        let hor_model_path = args
            .hor_model
            .clone()
            .unwrap_or_else(|| std::path::PathBuf::from(&cls_cfg.models.hor_score));
        let k_model_path = args
            .k_model
            .clone()
            .unwrap_or_else(|| std::path::PathBuf::from(&cls_cfg.models.k_pred));
        let hor_model = RandomForest::load_json(&hor_model_path)
            .with_context(|| format!("loading HOR-score model {:?}", hor_model_path))?;
        let k_model = RandomForest::load_json(&k_model_path).ok();
        if k_model.is_none() {
            warn!(
                "k-predictor model {:?} not loaded — k-recovery disabled",
                k_model_path
            );
        }
        let platt = cls_cfg.platt();

        let mut features: Vec<FeatureRow> = ok_records
            .par_iter()
            .zip(results.par_iter())
            .map(|(rec, kr)| build_features(rec, kr))
            .collect();

        if !args.no_homology {
            let mm_cfg = MonomerModelConfig::default();
            let mut probes: Vec<(usize, usize)> = Vec::with_capacity(features.len() * 2);
            for (i, f) in features.iter().enumerate() {
                if f.d1 > 0 {
                    probes.push((i, f.d1));
                }
                if f.family_founder_d > 0 && f.family_founder_d != f.d1 {
                    probes.push((i, f.family_founder_d));
                }
            }
            let probe_results: Vec<((usize, usize), Option<f64>)> = probes
                .par_iter()
                .map(|&(i, p)| {
                    let h = probe_period(ok_records[i], p, &mm_cfg).map(|(h, _, _)| h);
                    ((i, p), h)
                })
                .collect();
            for ((i, p), h) in probe_results {
                if let Some(h) = h {
                    if p == features[i].d1 {
                        features[i].h_d1 = h;
                    }
                    if p == features[i].family_founder_d {
                        features[i].h_founder = h;
                    }
                }
            }
            // h_founder defaults to h_d1 when no separate founder probe ran.
            for f in features.iter_mut() {
                if f.h_founder.is_nan() {
                    f.h_founder = f.h_d1;
                }
            }
        }

        let verdicts: Vec<(FeatureRow, Verdict)> = features
            .into_iter()
            .map(|mut f| {
                let v = run_classify(&mut f, &cls_cfg, &platt, &hor_model, k_model.as_ref());
                (f, v)
            })
            .collect();
        info!("classifier: applied to {} record(s)", verdicts.len());
        verdicts
    } else {
        Vec::new()
    };

    // Primary TSV.
    if let Some(parent) = args.out.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut out = std::fs::File::create(&args.out)
        .with_context(|| format!("creating {:?}", &args.out))?;
    let mut header = String::from(
        "case_id\tarray_length\tn_peaks_kept\
         \tmonomer_size\tscore\
         \tmonomer_size_2\tscore_2\
         \tmonomer_size_3\tscore_3",
    );
    if hor_enabled {
        header.push_str(
            "\thor_call\thor_founder\thor_multiplicity\thor_tile\
             \thor_family_size\thor_family_score\thor_jitter\thor_reason",
        );
    }
    if classify_enabled {
        header.push_str(
            "\thor_score\thor_score_raw\tverdict\
             \tfounder\tmultiplicity\ttile\tk_pred\trecovered\
             \th_d1\th_founder",
        );
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
        if hor_enabled {
            let hc = hor_classify(r, &hor_cfg);
            let f = hc
                .founder_bp
                .map(|x| x.to_string())
                .unwrap_or_else(|| na.into());
            let k = hc
                .multiplicity
                .map(|x| x.to_string())
                .unwrap_or_else(|| na.into());
            let t = hc.tile_bp.map(|x| x.to_string()).unwrap_or_else(|| na.into());
            line.push_str(&format!(
                "\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}",
                hc.verdict.as_str(),
                f,
                k,
                t,
                hc.family_size,
                hc.family_score,
                hc.jitter,
                hc.reason,
            ));
        }
        if classify_enabled {
            let (feat, verd) = &verdicts[idx];
            let fmt_opt = |o: &Option<usize>| {
                o.map(|x| x.to_string()).unwrap_or_else(|| na.into())
            };
            let fmt_h = |v: f64| {
                if v.is_nan() {
                    na.to_string()
                } else {
                    format!("{:.6}", v)
                }
            };
            line.push_str(&format!(
                "\t{:.10}\t{:.10}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                verd.hor_score,
                verd.hor_score_raw,
                verd.category.as_str(),
                fmt_opt(&verd.founder),
                fmt_opt(&verd.multiplicity),
                fmt_opt(&verd.tile),
                fmt_opt(&verd.k_pred),
                verd.recovered,
                fmt_h(feat.h_d1),
                fmt_h(feat.h_founder),
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
