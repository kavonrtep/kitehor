#!/usr/bin/env Rscript
# RF with the +2 homology features (h_d1, h_founder).

suppressPackageStartupMessages({ library(ranger) })

train_path <- "eval/training_data/master_features_train_h.tsv"
test_path  <- "eval/training_data/master_features_test_h.tsv"
outdir     <- "eval/reports/hor_model_rf_h"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

train <- read.table(train_path, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, na.strings = "NA")
test  <- read.table(test_path,  header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, na.strings = "NA")

DETECTABLE_THR <- 0.10
train$y <- factor(as.integer(train$hor_signal >= DETECTABLE_THR), levels = c(0, 1))
test$y  <- factor(as.integer(test$hor_signal  >= DETECTABLE_THR), levels = c(0, 1))

# Impute homology NAs with the column's training median (cleanest for RF).
for (col in c("h_d1", "h_founder")) {
  med <- median(train[[col]], na.rm = TRUE)
  train[[col]][is.na(train[[col]])] <- med
  test[[col]][is.na(test[[col]])]   <- med
  cat(sprintf("%s: %d NAs filled with median=%.4f (train)\n",
              col, sum(is.na(c(train[[col]], test[[col]]))), med))
}

feature_cols <- c(
  "s1", "s2", "s3", "s2_over_s1", "s3_over_s1",
  "family_size_best", "tile_founder_ratio", "tile_jitter",
  "d1", "log_d1_over_L",
  "distinct_kmers_per_bp", "kmer_entropy", "singletons_ratio",
  "h_d1", "h_founder"
)
cat(sprintf("\ntrain: %d (%d positive)\n", nrow(train), sum(train$y == 1)))
cat(sprintf("test:  %d (%d positive)\n\n", nrow(test),  sum(test$y == 1)))

set.seed(42)
fit <- ranger(
  formula = as.formula(paste("y ~", paste(feature_cols, collapse = " + "))),
  data = train, num.trees = 500, importance = "permutation",
  probability = TRUE, verbose = FALSE
)
cat(sprintf("OOB Brier: %.4f\n", fit$prediction.error))

pred <- predict(fit, data = test)
test$score <- pred$predictions[, "1"]

imp <- data.frame(feature = names(fit$variable.importance),
                  importance = unname(fit$variable.importance))
imp <- imp[order(-imp$importance), ]
write.table(imp, file.path(outdir, "feature_importance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n=== Feature importance ===\n"); print(imp); cat("\n")

# PR/ROC
get_metrics <- function(score, label, thresholds) {
  out <- data.frame(threshold = thresholds, tp = NA_integer_,
                    fp = NA_integer_, tn = NA_integer_, fn = NA_integer_,
                    precision = NA_real_, recall = NA_real_,
                    f1 = NA_real_, fpr = NA_real_)
  lab <- as.integer(as.character(label))
  for (i in seq_along(thresholds)) {
    pr <- score >= thresholds[i]
    tp <- sum(pr & lab == 1); fp <- sum(pr & lab == 0)
    tn <- sum(!pr & lab == 0); fn <- sum(!pr & lab == 1)
    out$tp[i] <- tp; out$fp[i] <- fp; out$tn[i] <- tn; out$fn[i] <- fn
    out$precision[i] <- if (tp+fp==0) NA else tp/(tp+fp)
    out$recall[i]    <- if (tp+fn==0) NA else tp/(tp+fn)
    out$f1[i] <- if (is.na(out$precision[i])||is.na(out$recall[i])
                    ||out$precision[i]+out$recall[i]==0) NA
                 else 2*out$precision[i]*out$recall[i]/
                      (out$precision[i]+out$recall[i])
    out$fpr[i] <- if (fp+tn==0) NA else fp/(fp+tn)
  }
  out
}
pr <- get_metrics(test$score, test$y, seq(0, 1, 0.01))
write.table(pr, file.path(outdir, "pr_curve_overall.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

auc <- function(x, y) {
  ord <- order(x); x <- x[ord]; y <- y[ord]
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(NA)
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}
pr_auc  <- auc(rev(pr$recall), rev(pr$precision))
roc_auc <- auc(pr$fpr, pr$recall)
cat(sprintf("PR  AUC (held-out): %.4f\n", pr_auc))
cat(sprintf("ROC AUC (held-out): %.4f\n\n", roc_auc))

best_f1_idx <- which.max(pr$f1)
hi_idx <- which(pr$precision >= 0.95)[1]
lo_idx <- which(pr$precision >= 0.50)[1]
chosen <- data.frame(
  name = c("t_low (precision≥0.50)",
           "t_high (precision≥0.95)", "t_max_F1"),
  threshold = c(pr$threshold[lo_idx], pr$threshold[hi_idx],
                pr$threshold[best_f1_idx]),
  precision = c(pr$precision[lo_idx], pr$precision[hi_idx],
                pr$precision[best_f1_idx]),
  recall = c(pr$recall[lo_idx], pr$recall[hi_idx],
             pr$recall[best_f1_idx]),
  f1 = c(pr$f1[lo_idx], pr$f1[hi_idx], pr$f1[best_f1_idx])
)
write.table(chosen, file.path(outdir, "thresholds.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("=== Chosen thresholds ===\n"); print(chosen); cat("\n")

# Per-stratum
chosen_thr <- chosen$threshold[chosen$name == "t_max_F1"]
test$pred <- as.integer(test$score >= chosen_thr)
test$y_int <- as.integer(as.character(test$y))
per_stratum <- data.frame(
  stratum = sort(unique(test$stratum)),
  n = 0L, n_pos = 0L,
  tp = 0L, fp = 0L, tn = 0L, fn = 0L,
  precision = NA_real_, recall = NA_real_, f1 = NA_real_
)
for (i in seq_along(per_stratum$stratum)) {
  sub <- test[test$stratum == per_stratum$stratum[i], ]
  per_stratum$n[i] <- nrow(sub)
  per_stratum$n_pos[i] <- sum(sub$y_int == 1)
  tp <- sum(sub$pred == 1 & sub$y_int == 1)
  fp <- sum(sub$pred == 1 & sub$y_int == 0)
  tn <- sum(sub$pred == 0 & sub$y_int == 0)
  fn <- sum(sub$pred == 0 & sub$y_int == 1)
  per_stratum$tp[i] <- tp; per_stratum$fp[i] <- fp
  per_stratum$tn[i] <- tn; per_stratum$fn[i] <- fn
  per_stratum$precision[i] <- if (tp+fp==0) NA else tp/(tp+fp)
  per_stratum$recall[i]    <- if (tp+fn==0) NA else tp/(tp+fn)
  per_stratum$f1[i] <- if (is.na(per_stratum$precision[i])||is.na(per_stratum$recall[i])
                          ||per_stratum$precision[i]+per_stratum$recall[i]==0) NA
                       else 2*per_stratum$precision[i]*per_stratum$recall[i]/
                            (per_stratum$precision[i]+per_stratum$recall[i])
}
write.table(per_stratum, file.path(outdir, "per_stratum.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("=== Per-stratum at threshold=%.2f ===\n", chosen_thr))
print(per_stratum); cat("\n")

saveRDS(fit, file.path(outdir, "ranger_model.rds"))

sink(file.path(outdir, "summary.txt"))
cat(sprintf("RF + homology — binary on hor_signal >= %.2f\n", DETECTABLE_THR))
cat(sprintf("Features: %d\n", length(feature_cols)))
cat(sprintf("Train: %d (%d positive) | Test: %d (%d positive)\n",
            nrow(train), sum(train$y == 1), nrow(test), sum(test$y == 1)))
cat(sprintf("PR  AUC: %.4f  | ROC AUC: %.4f\n\n", pr_auc, roc_auc))
cat("Importance:\n"); print(imp)
cat("\nThresholds:\n"); print(chosen)
cat("\nPer-stratum at t_max_F1:\n"); print(per_stratum)
sink()
cat(sprintf("\nOutputs in %s\n", outdir))
