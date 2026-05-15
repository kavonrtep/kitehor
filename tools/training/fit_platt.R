#!/usr/bin/env Rscript
# Platt scaling: fit a logistic regression that maps the RF score
# (uncalibrated probability) to a calibrated posterior probability.
#
# We use the RF's OUT-OF-BAG predictions on the training set as
# "honest" scores (each record's prediction comes from trees that
# didn't see it). This avoids overfitting that would result from
# calibrating on in-sample scores.
#
# Output:
#   eval/reports/hor_model_rf_h/platt_calibration.rds
#   eval/reports/hor_model_rf_h/calibration_diagnostics.tsv
#   eval/reports/hor_model_rf_h/thresholds_calibrated.tsv

suppressPackageStartupMessages({ library(ranger) })

fit <- readRDS("eval/reports/hor_model_rf_h/ranger_model.rds")
train <- read.table("eval/training_data/master_features_train_h.tsv",
                    header = TRUE, sep = "\t", na.strings = "NA",
                    stringsAsFactors = FALSE)

# Impute homology medians (same as training).
med_h_d1 <- median(train$h_d1, na.rm = TRUE)
med_h_founder <- median(train$h_founder, na.rm = TRUE)
train$h_d1[is.na(train$h_d1)] <- med_h_d1
train$h_founder[is.na(train$h_founder)] <- med_h_founder

DETECTABLE_THR <- 0.10
y_train <- as.integer(train$hor_signal >= DETECTABLE_THR)

# OOB scores from ranger (probability mode). Each training record was
# predicted by trees that didn't include it as a bootstrap sample.
oob_scores <- fit$predictions[, "1"]

# Clip to avoid logit infinity at edges.
EPS <- 1e-4
oob_clip <- pmin(pmax(oob_scores, EPS), 1 - EPS)
logit_score <- log(oob_clip / (1 - oob_clip))

# Platt scaling: P(y=1 | s) = sigmoid(a + b·logit(s)).
platt <- glm(y_train ~ logit_score, family = binomial(link = "logit"))
cat("=== Platt scaling fit ===\n")
print(summary(platt))
cat("\n")

# Save fit + medians (needed at inference time for imputation).
saveRDS(list(platt = platt,
             med_h_d1 = med_h_d1,
             med_h_founder = med_h_founder),
        "eval/reports/hor_model_rf_h/platt_calibration.rds")
cat("Saved calibration → eval/reports/hor_model_rf_h/platt_calibration.rds\n\n")

# Apply Platt to OOB scores and bin to check calibration.
calibrated_oob <- predict(platt, newdata = data.frame(logit_score = logit_score),
                          type = "response")
bins <- cut(calibrated_oob, breaks = seq(0, 1, by = 0.1),
            include.lowest = TRUE)
diag <- aggregate(cbind(
    n = rep(1, length(calibrated_oob)),
    cal = calibrated_oob,
    raw = oob_scores,
    actual = y_train
  ) ~ bins, FUN = mean)
diag$n <- as.numeric(table(bins))[match(diag$bins, names(table(bins)))]
write.table(diag, "eval/reports/hor_model_rf_h/calibration_diagnostics.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("=== Calibration diagnostic (OOB training) ===\n")
cat("Bin            n   raw_mean  cal_mean  actual_rate\n")
for (i in seq_len(nrow(diag))) {
  cat(sprintf("%-12s %5d  %.4f    %.4f    %.4f\n",
              diag$bins[i], diag$n[i], diag$raw[i],
              diag$cal[i], diag$actual[i]))
}
cat("\n")

# Re-tune thresholds on the held-out original test set
# (seeds 901+902, never used for RF training).
test <- read.table("eval/training_data/master_features_test_h.tsv",
                   header = TRUE, sep = "\t", na.strings = "NA",
                   stringsAsFactors = FALSE)
test$h_d1[is.na(test$h_d1)] <- med_h_d1
test$h_founder[is.na(test$h_founder)] <- med_h_founder
test$y <- as.integer(test$hor_signal >= DETECTABLE_THR)
test_pred <- predict(fit, data = test)$predictions[, "1"]
test_pred_clip <- pmin(pmax(test_pred, EPS), 1 - EPS)
test$score_raw <- test_pred
test$score <- predict(platt,
    newdata = data.frame(logit_score = log(test_pred_clip / (1 - test_pred_clip))),
    type = "response")

# PR curve.
get_metrics <- function(score, y, thresholds) {
  res <- data.frame(threshold = thresholds, precision = NA_real_,
                    recall = NA_real_, f1 = NA_real_)
  for (i in seq_along(thresholds)) {
    p <- score >= thresholds[i]
    tp <- sum(p & y == 1); fp <- sum(p & y == 0); fn <- sum(!p & y == 1)
    res$precision[i] <- if (tp+fp==0) NA else tp/(tp+fp)
    res$recall[i]    <- if (tp+fn==0) NA else tp/(tp+fn)
    res$f1[i] <- if (is.na(res$precision[i])||is.na(res$recall[i])
                    ||(res$precision[i]+res$recall[i])==0) NA
                 else 2*res$precision[i]*res$recall[i]/
                      (res$precision[i]+res$recall[i])
  }
  res
}
pr <- get_metrics(test$score, test$y, seq(0, 1, by = 0.01))
hi_idx  <- which(pr$precision >= 0.95)[1]
lo_idx  <- which(pr$precision >= 0.50)[1]
f1_idx  <- which.max(pr$f1)

chosen <- data.frame(
  name = c("t_low_calibrated (precision≥0.50)",
           "t_high_calibrated (precision≥0.95)",
           "t_max_F1_calibrated"),
  threshold = c(pr$threshold[lo_idx],
                pr$threshold[hi_idx],
                pr$threshold[f1_idx]),
  precision = c(pr$precision[lo_idx],
                pr$precision[hi_idx],
                pr$precision[f1_idx]),
  recall = c(pr$recall[lo_idx],
             pr$recall[hi_idx],
             pr$recall[f1_idx]),
  f1 = c(pr$f1[lo_idx], pr$f1[hi_idx], pr$f1[f1_idx])
)
write.table(chosen, "eval/reports/hor_model_rf_h/thresholds_calibrated.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("=== Calibrated thresholds on held-out test ===\n")
print(chosen)
