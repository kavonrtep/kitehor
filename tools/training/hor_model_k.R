#!/usr/bin/env Rscript
# Train a k-predictor: RF regression on truth.hor_order, restricted to
# truth-HOR records (k >= 2). Same 15-feature set as the hor_score model.
#
# The motivation: when the hor_score model says HOR but kite didn't find
# a clean integer-multiple family, we have no (founder, tile, k). A
# direct k-prediction can fill this gap — knowing predicted k, we can
# then look for d1·k or d1/k in the kite peaks.
#
# Output:
#   eval/reports/hor_model_k/
#     ranger_model_k.rds              — trained regressor
#     feature_importance.tsv
#     pr_curve_overall.tsv           — N/A here, classification-style metrics
#     accuracy_by_truth_k.tsv         — per-truth-k accuracy
#     predictions_cv.tsv              — held-out predictions
#     summary.txt

suppressPackageStartupMessages({ library(ranger) })

train_path <- "eval/training_data/master_features_train_h.tsv"
test_path  <- "eval/training_data/master_features_test_h.tsv"
outdir     <- "eval/reports/hor_model_k"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

train <- read.table(train_path, header = TRUE, sep = "\t",
                    na.strings = "NA", stringsAsFactors = FALSE)
test  <- read.table(test_path,  header = TRUE, sep = "\t",
                    na.strings = "NA", stringsAsFactors = FALSE)

# Pull hor_order from the original truth file (it's not yet in features).
get_hor_order <- function(features_dir_root) {
  seeds <- list.files(features_dir_root, pattern = "^sim_seed", full.names = TRUE)
  all <- list()
  for (sd in seeds) {
    tp <- file.path(sd, "truth.tsv")
    if (!file.exists(tp)) next
    t <- read.table(tp, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE)
    all[[sd]] <- t[, c("case_id", "hor_order")]
  }
  do.call(rbind, all)
}

ho_train <- do.call(rbind, lapply(c(101, 102, 103, 104, 105), function(s) {
  read.table(sprintf("eval/training_data/sim_seed%d/truth.tsv", s),
             header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, c("case_id", "hor_order")]
}))
ho_test <- do.call(rbind, lapply(c(901, 902), function(s) {
  read.table(sprintf("eval/training_data/sim_seed%d/truth.tsv", s),
             header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, c("case_id", "hor_order")]
}))
# Some case_ids may appear in multiple seeds; we need to match by (seed,case_id).
# But seeds are concatenated and case_id is unique-ish per file; the master TSV
# preserves order so just take the column-wise hor_order per row.
# Simpler: read truth per seed and concatenate IN THE SAME ORDER as features.
truth_train <- do.call(rbind, lapply(c(101, 102, 103, 104, 105), function(s) {
  read.table(sprintf("eval/training_data/sim_seed%d/truth.tsv", s),
             header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, c("case_id", "hor_order")]
}))
truth_test <- do.call(rbind, lapply(c(901, 902), function(s) {
  read.table(sprintf("eval/training_data/sim_seed%d/truth.tsv", s),
             header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, c("case_id", "hor_order")]
}))
stopifnot(nrow(truth_train) == nrow(train))
stopifnot(nrow(truth_test) == nrow(test))
train$k_truth <- truth_train$hor_order
test$k_truth  <- truth_test$hor_order

# Impute h_* NAs (use training medians).
med_h_d1 <- median(train$h_d1, na.rm = TRUE)
med_h_founder <- median(train$h_founder, na.rm = TRUE)
train$h_d1[is.na(train$h_d1)] <- med_h_d1
train$h_founder[is.na(train$h_founder)] <- med_h_founder
test$h_d1[is.na(test$h_d1)] <- med_h_d1
test$h_founder[is.na(test$h_founder)] <- med_h_founder

# Filter to truth-HOR (k >= 2) for training.
train_hor <- train[train$k_truth >= 2, ]
test_hor  <- test[test$k_truth >= 2, ]
cat(sprintf("train (HOR-only): %d rows, k range = [%d, %d]\n",
            nrow(train_hor), min(train_hor$k_truth), max(train_hor$k_truth)))
cat(sprintf("test  (HOR-only): %d rows, k range = [%d, %d]\n\n",
            nrow(test_hor),  min(test_hor$k_truth),  max(test_hor$k_truth)))

cat("k distribution in train:\n"); print(table(train_hor$k_truth))
cat("\nk distribution in test:\n");  print(table(test_hor$k_truth)); cat("\n")

feature_cols <- c(
  "s1", "s2", "s3", "s2_over_s1", "s3_over_s1",
  "family_size_best", "tile_founder_ratio", "tile_jitter",
  "d1", "d2", "d3", "log_d1_over_L",
  "d2_over_d1", "d3_over_d1", "max_d_top3_over_min_d_top3",
  "distinct_kmers_per_bp", "kmer_entropy", "singletons_ratio",
  "h_d1", "h_founder"
)
set.seed(43)
fit_k <- ranger(
  formula = as.formula(paste("k_truth ~", paste(feature_cols, collapse = " + "))),
  data = train_hor, num.trees = 500,
  importance = "permutation", verbose = FALSE
)
cat(sprintf("OOB MSE (training): %.4f  (RMSE %.3f)\n",
            fit_k$prediction.error, sqrt(fit_k$prediction.error)))

# Held-out predictions.
test_hor$k_pred_cont <- predict(fit_k, data = test_hor)$predictions
test_hor$k_pred <- round(test_hor$k_pred_cont)
test_hor$k_pred <- pmax(test_hor$k_pred, 2)  # enforce k >= 2 in HOR-only setting

# Overall accuracy.
acc <- mean(test_hor$k_pred == test_hor$k_truth)
acc_within1 <- mean(abs(test_hor$k_pred - test_hor$k_truth) <= 1)
cat(sprintf("\nHeld-out accuracy (k_pred == k_truth):       %.3f\n", acc))
cat(sprintf("Held-out accuracy (|k_pred - k_truth| <= 1):  %.3f\n", acc_within1))
cat(sprintf("Held-out RMSE on continuous k:                %.3f\n",
            sqrt(mean((test_hor$k_pred_cont - test_hor$k_truth)^2))))

cat("\n=== Confusion (truth k vs predicted k) ===\n")
print(table(truth_k = test_hor$k_truth, pred_k = test_hor$k_pred))

# Per-truth-k accuracy
per_k <- aggregate(k_pred == k_truth ~ k_truth, data = test_hor,
                   FUN = function(x) c(n = length(x),
                                        accuracy = round(mean(x), 3),
                                        accuracy_within1 = round(mean(abs(x) < 2), 3)))
cat("\n=== Per-truth-k accuracy ===\n")
print(per_k)

# Per-stratum predictions
per_stratum <- aggregate(cbind(k_truth, k_pred_cont, k_pred) ~ stratum,
                          data = test_hor,
                          FUN = function(x) round(c(n = length(x),
                                                     mean = mean(x),
                                                     median = median(x),
                                                     min = min(x),
                                                     max = max(x)), 2))
cat("\n=== Per-stratum k summary ===\n")
print(per_stratum)

# Feature importance
imp <- data.frame(feature = names(fit_k$variable.importance),
                  importance = unname(fit_k$variable.importance))
imp <- imp[order(-imp$importance), ]
write.table(imp, file.path(outdir, "feature_importance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n=== Feature importance ===\n"); print(imp); cat("\n")

# Save model + predictions
saveRDS(fit_k, file.path(outdir, "ranger_model_k.rds"))
out_pred <- test_hor[, c("case_id", "stratum", "k_truth",
                          "k_pred_cont", "k_pred")]
write.table(out_pred, file.path(outdir, "predictions_holdout.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sink(file.path(outdir, "summary.txt"))
cat(sprintf("k-predictor (RF regression on truth.hor_order, k >= 2)\n"))
cat(sprintf("Train: %d rows, Test: %d rows\n", nrow(train_hor), nrow(test_hor)))
cat(sprintf("OOB MSE: %.3f (RMSE %.3f)\n",
            fit_k$prediction.error, sqrt(fit_k$prediction.error)))
cat(sprintf("Held-out accuracy: %.3f exact, %.3f within ±1\n",
            acc, acc_within1))
cat(sprintf("Held-out RMSE on continuous k: %.3f\n",
            sqrt(mean((test_hor$k_pred_cont - test_hor$k_truth)^2))))
cat("\nImportance:\n"); print(imp)
cat("\nPer-truth-k accuracy:\n"); print(per_k)
sink()
cat(sprintf("\nOutputs in %s\n", outdir))
