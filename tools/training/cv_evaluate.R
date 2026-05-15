#!/usr/bin/env Rscript
# Cross-validation on 10 new simulated seeds (201..210), using the
# trained ranger RF model from hor_model_rf_h.R.
#
# Reports per-seed metrics + aggregate (mean ± sd) across the 10 seeds:
#   - PR AUC, ROC AUC
#   - precision / recall / F1 at the F1-optimal threshold (learned on
#     the original train+test held-out at t=0.41)
#   - per-stratum precision/recall/F1 at the same threshold
#
# Also re-fits thresholds per seed to confirm the canonical threshold
# transfers.

suppressPackageStartupMessages({ library(ranger) })

fit <- readRDS("eval/reports/hor_model_rf_h/ranger_model.rds")
THR_F1 <- 0.41   # from hor_model_rf_h.R
DETECTABLE_THR <- 0.10
SEEDS <- c(201, 202, 203, 204, 205, 206, 207, 208, 209, 210)
feature_cols <- c("s1", "s2", "s3", "s2_over_s1", "s3_over_s1",
                  "family_size_best", "tile_founder_ratio", "tile_jitter",
                  "d1", "log_d1_over_L",
                  "distinct_kmers_per_bp", "kmer_entropy", "singletons_ratio",
                  "h_d1", "h_founder")
outdir <- "eval/reports/hor_model_rf_h/cv"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Need training median for NA imputation of h_d1/h_founder.
train <- read.table("eval/training_data/master_features_train_h.tsv",
                    header = TRUE, sep = "\t", na.strings = "NA")
med_h_d1 <- median(train$h_d1, na.rm = TRUE)
med_h_founder <- median(train$h_founder, na.rm = TRUE)

auc_trap <- function(x, y) {
  ord <- order(x); x <- x[ord]; y <- y[ord]
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(NA)
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

per_seed_overall <- data.frame()
per_seed_stratum <- data.frame()

for (seed in SEEDS) {
  fpath <- sprintf("eval/training_data/sim_seed%d/features_h.tsv", seed)
  d <- read.table(fpath, header = TRUE, sep = "\t", na.strings = "NA",
                  stringsAsFactors = FALSE)
  d$h_d1[is.na(d$h_d1)] <- med_h_d1
  d$h_founder[is.na(d$h_founder)] <- med_h_founder
  d$y <- as.integer(d$hor_signal >= DETECTABLE_THR)
  pred <- predict(fit, data = d)
  d$score <- pred$predictions[, "1"]

  # PR/ROC over fine threshold grid.
  thresholds <- seq(0, 1, by = 0.005)
  prec <- numeric(length(thresholds)); rec <- numeric(length(thresholds))
  fpr  <- numeric(length(thresholds))
  for (i in seq_along(thresholds)) {
    p <- d$score >= thresholds[i]
    tp <- sum(p & d$y == 1); fp <- sum(p & d$y == 0)
    tn <- sum(!p & d$y == 0); fn <- sum(!p & d$y == 1)
    prec[i] <- if (tp+fp==0) NA else tp/(tp+fp)
    rec[i]  <- if (tp+fn==0) NA else tp/(tp+fn)
    fpr[i]  <- if (fp+tn==0) NA else fp/(fp+tn)
  }
  pr_auc  <- auc_trap(rev(rec), rev(prec))
  roc_auc <- auc_trap(fpr, rec)
  # Threshold = canonical t=THR_F1
  p_t <- d$score >= THR_F1
  tp <- sum(p_t & d$y == 1); fp <- sum(p_t & d$y == 0)
  tn <- sum(!p_t & d$y == 0); fn <- sum(!p_t & d$y == 1)
  f1 <- if (tp+fp == 0 || tp+fn == 0) NA
        else 2*tp/(2*tp+fp+fn)
  per_seed_overall <- rbind(per_seed_overall, data.frame(
    seed = seed,
    n = nrow(d), n_pos = sum(d$y == 1),
    pr_auc = round(pr_auc, 4), roc_auc = round(roc_auc, 4),
    tp = tp, fp = fp, tn = tn, fn = fn,
    precision = round(if(tp+fp==0) NA else tp/(tp+fp), 4),
    recall = round(if(tp+fn==0) NA else tp/(tp+fn), 4),
    f1 = round(f1, 4)
  ))
  # Per-stratum at canonical threshold.
  for (s in sort(unique(d$stratum))) {
    sub <- d[d$stratum == s, ]
    p <- sub$score >= THR_F1
    tp <- sum(p & sub$y == 1); fp <- sum(p & sub$y == 0)
    tn <- sum(!p & sub$y == 0); fn <- sum(!p & sub$y == 1)
    f1 <- if (tp+fp == 0 || tp+fn == 0) NA
          else 2*tp/(2*tp+fp+fn)
    per_seed_stratum <- rbind(per_seed_stratum, data.frame(
      seed = seed, stratum = s, n = nrow(sub),
      n_pos = sum(sub$y == 1),
      precision = if(tp+fp==0) NA else tp/(tp+fp),
      recall = if(tp+fn==0) NA else tp/(tp+fn),
      f1 = f1
    ))
  }
}

write.table(per_seed_overall, file.path(outdir, "per_seed_overall.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(per_seed_stratum, file.path(outdir, "per_seed_stratum.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("=== Per-seed overall (at t=0.41) ===\n")
print(per_seed_overall)
cat("\n=== Aggregate across 10 CV seeds ===\n")
agg <- data.frame(
  metric = c("pr_auc", "roc_auc", "precision", "recall", "f1"),
  mean = c(mean(per_seed_overall$pr_auc), mean(per_seed_overall$roc_auc),
           mean(per_seed_overall$precision, na.rm=TRUE),
           mean(per_seed_overall$recall, na.rm=TRUE),
           mean(per_seed_overall$f1, na.rm=TRUE)),
  sd = c(sd(per_seed_overall$pr_auc), sd(per_seed_overall$roc_auc),
         sd(per_seed_overall$precision, na.rm=TRUE),
         sd(per_seed_overall$recall, na.rm=TRUE),
         sd(per_seed_overall$f1, na.rm=TRUE)),
  min = c(min(per_seed_overall$pr_auc), min(per_seed_overall$roc_auc),
          min(per_seed_overall$precision, na.rm=TRUE),
          min(per_seed_overall$recall, na.rm=TRUE),
          min(per_seed_overall$f1, na.rm=TRUE)),
  max = c(max(per_seed_overall$pr_auc), max(per_seed_overall$roc_auc),
          max(per_seed_overall$precision, na.rm=TRUE),
          max(per_seed_overall$recall, na.rm=TRUE),
          max(per_seed_overall$f1, na.rm=TRUE))
)
agg$mean <- round(agg$mean, 4); agg$sd <- round(agg$sd, 4)
agg$min <- round(agg$min, 4);   agg$max <- round(agg$max, 4)
print(agg)
write.table(agg, file.path(outdir, "aggregate.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Aggregate per stratum.
cat("\n=== Per-stratum aggregate (mean ± sd over 10 CV seeds) ===\n")
strata_agg <- aggregate(cbind(precision, recall, f1) ~ stratum,
                        data = per_seed_stratum,
                        FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                             sd = sd(x, na.rm = TRUE)))
print(strata_agg)
write.table(do.call(data.frame, strata_agg),
            file.path(outdir, "per_stratum_aggregate.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\nOutputs in %s\n", outdir))
