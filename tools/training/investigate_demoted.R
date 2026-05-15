#!/usr/bin/env Rscript
# Investigate the residual demoted cases — those that were predicted as
# HOR by the RF (score >= 0.60) but had no kite family (family_founder_d=0)
# AND were NOT recovered by the k-predictor.
#
# Stratify by:
#   - stratum
#   - truth.monomer_len (founder length)
#   - truth.hor_order   (k)
#   - truth.hor_signal
#   - features (kite scores, diversity, homology)
#
# Look for systematic patterns: is there a "type" of HOR the pipeline
# consistently fails to call?

# Read predictions + join truth metadata.
pred <- read.table("eval/reports/hor_model_rf_h/verdict/cv_pool_predictions_v3.tsv",
                   header = TRUE, sep = "\t", na.strings = "NA",
                   stringsAsFactors = FALSE)
# Cross-reference truth metadata by reading truth.tsv per seed and matching.
seeds <- c(201, 202, 203, 204, 205, 206, 207, 208, 209, 210)
truth_all <- do.call(rbind, lapply(seeds, function(s) {
  t <- read.table(sprintf("eval/training_data/sim_seed%d/truth.tsv", s),
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  t$seed <- s
  t
}))
# pred has case_ids that repeat across seeds (because params.tsv is the same).
# Matching requires (seed, case_id). The pred row order matches the
# concatenation order of features tsvs. Read seed-by-seed counts.
pred$seed <- rep(seeds, times = sapply(seeds, function(s) {
  nrow(read.table(sprintf("eval/training_data/sim_seed%d/features_h.tsv", s),
                  header = TRUE, sep = "\t"))
}))
# Drop hor_signal from pred (it's already in truth_all) to avoid merge ambiguity.
pred$hor_signal <- NULL
pred$hor_binary <- NULL
# Keep needed truth columns only
truth_keep <- truth_all[, c("seed", "case_id", "monomer_len", "hor_order",
                             "hor_signal", "n_blocks", "array_length",
                             "submono_k")]
colnames(truth_keep)[colnames(truth_keep) == "array_length"] <- "truth_array_length"
df <- merge(pred, truth_keep, by = c("seed", "case_id"))
cat(sprintf("merged: %d rows\n", nrow(df)))

# Mark each case's bucket relative to the verdict pipeline.
df$truth_hor <- as.integer(df$hor_order > 1)
df$truth_detectable <- as.integer(df$hor_signal >= 0.10)

# True HOR cases that *should* have been called HOR (truth.hor_signal >= 0.10)
# but landed in unresolved (i.e., the RF said HOR but kite found no family
# and k-recovery failed).
demoted <- df[df$hor_score >= 0.60 &
              (df$family_founder_d == 0 | df$family_founder_d == df$family_tile_d) &
              !df$recovered, ]
cat(sprintf("residual demoted (not recovered): %d rows\n", nrow(demoted)))
cat("Of these, true-HOR (hor_signal >= 0.10):",
    sum(demoted$hor_signal >= 0.10), "/", nrow(demoted), "\n\n")

# === 1. By stratum ===
cat("=== Demoted cases by stratum ===\n")
strat <- aggregate(cbind(n = rep(1, nrow(demoted)),
                          true_hor = demoted$hor_signal >= 0.10,
                          true_not_hor = demoted$hor_signal < 0.10) ~ stratum,
                    data = demoted, FUN = sum)
strat$pct_true_hor <- round(100 * strat$true_hor / strat$n, 1)
print(strat)
cat("\n")

# === 2. By monomer_len bucket ===
cat("=== Demoted cases by truth.monomer_len ===\n")
demoted$m_bucket <- cut(demoted$monomer_len,
                        breaks = c(0, 50, 100, 200, 500, 1000, 3000),
                        right = TRUE)
mb <- aggregate(cbind(n = rep(1, nrow(demoted)),
                       true_hor = demoted$hor_signal >= 0.10) ~ m_bucket,
                 data = demoted, FUN = sum)
mb$pct_true_hor <- round(100 * mb$true_hor / mb$n, 1)
print(mb)
cat("\n")

# === 3. By truth.hor_order ===
cat("=== Demoted cases by truth.hor_order ===\n")
ho <- aggregate(cbind(n = rep(1, nrow(demoted)),
                       true_hor = demoted$hor_signal >= 0.10) ~ hor_order,
                 data = demoted, FUN = sum)
ho$pct_true_hor <- round(100 * ho$true_hor / ho$n, 1)
print(ho)
cat("\n")

# === 4. By hor_signal bucket ===
cat("=== Demoted cases by truth.hor_signal ===\n")
demoted$hs_bucket <- cut(demoted$hor_signal,
                          breaks = c(-0.1, 0.0001, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50),
                          labels = c("0", "(0,0.05]", "(0.05,0.10]",
                                     "(0.10,0.15]", "(0.15,0.20]",
                                     "(0.20,0.30]", ">0.30"),
                          include.lowest = TRUE)
hs <- as.data.frame(table(demoted$hs_bucket))
colnames(hs) <- c("hs_bucket", "n")
print(hs)
cat("\n")

# === 5. Stratum × hor_order cross-tab ===
cat("=== Demoted: stratum × hor_order ===\n")
print(addmargins(table(stratum = demoted$stratum, hor_order = demoted$hor_order)))
cat("\n")

# === 6. Feature distributions ===
cat("=== Feature distributions of demoted cases (true HOR vs not) ===\n")
true_hor <- demoted[demoted$hor_signal >= 0.10, ]
true_not <- demoted[demoted$hor_signal < 0.10, ]
feats_to_summarize <- c("hor_score", "k_pred", "s1", "s2", "s3",
                        "family_size_best", "tile_jitter",
                        "d1", "d2", "d3",
                        "max_d_top3_over_min_d_top3",
                        "h_d1", "h_founder",
                        "kmer_entropy")
for (f in feats_to_summarize) {
  th <- true_hor[[f]]; tn <- true_not[[f]]
  cat(sprintf("%-30s  true_HOR median=%-8.3f  true_not median=%-8.3f\n",
              f, median(th, na.rm=TRUE), median(tn, na.rm=TRUE)))
}
cat("\n")

# === 7. Compare RECOVERED cases to demoted: what discriminated them? ===
cat("=== Recovered vs demoted (when family was missing): which features differ? ===\n")
recovered <- df[df$recovered, ]
also_demoted <- df[df$hor_score >= 0.60 &
                   (df$family_founder_d == 0 | df$family_founder_d == df$family_tile_d) &
                   !df$recovered, ]
cat(sprintf("Recovered: %d, Still-demoted: %d\n", nrow(recovered), nrow(also_demoted)))
for (f in c("hor_score", "k_pred", "d1", "d2", "d3",
            "max_d_top3_over_min_d_top3", "s2_over_s1", "h_d1")) {
  cat(sprintf("%-30s  recovered=%-8.3f  demoted=%-8.3f\n",
              f,
              median(recovered[[f]], na.rm = TRUE),
              median(also_demoted[[f]], na.rm = TRUE)))
}
cat("\n")

# Save the demoted records.
write.table(demoted, "eval/reports/hor_model_rf_h/verdict/demoted_inspect.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("saved → eval/reports/hor_model_rf_h/verdict/demoted_inspect.tsv\n")
