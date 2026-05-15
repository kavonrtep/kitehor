#!/usr/bin/env Rscript
# Apply the trained ranger RF model + 4-category verdict logic.
#
# Verdict policy (in order):
#   1. kite found no peaks (s1 == 0 OR d1 == 0)         → no_signal
#   2. score >= t_high (default 0.60)                    → hor
#   3. score <  t_low  (default 0.20)                    → tandem
#   4. otherwise                                          → unresolved
#
# Output: input features + score + verdict + (founder, multiplicity, tile)
# when applicable, in a per-record TSV.
#
# Usage:
#   Rscript predict_verdict.R --features <in.tsv> --out <out.tsv>
#       [--model <path>]   default: eval/reports/hor_model_rf_h/ranger_model.rds
#       [--t-low 0.20] [--t-high 0.60]

suppressPackageStartupMessages({ library(ranger) })

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(key, default = NULL) {
  i <- which(args == key)
  if (length(i) == 0) return(default)
  args[i + 1]
}
features_path <- parse_arg("--features")
out_path      <- parse_arg("--out")
model_path    <- parse_arg("--model",
                           "eval/reports/hor_model_rf_h/ranger_model.rds")
train_path    <- parse_arg("--train",
                           "eval/training_data/master_features_train_h.tsv")
T_LOW         <- as.numeric(parse_arg("--t-low",  "0.15"))
T_HIGH        <- as.numeric(parse_arg("--t-high", "0.71"))
CALIBRATION   <- parse_arg("--calibration",
                           "eval/reports/hor_model_rf_h/platt_calibration.rds")
if (is.null(features_path) || is.null(out_path)) {
  stop("usage: predict_verdict.R --features <in.tsv> --out <out.tsv>")
}

# Read training data once to compute h_* imputation medians (the RF was
# trained with NA-imputed homology values; we apply the same imputation
# at inference time so the model sees the same feature distribution).
train <- read.table(train_path, header = TRUE, sep = "\t",
                    na.strings = "NA", stringsAsFactors = FALSE)
med_h_d1 <- median(train$h_d1, na.rm = TRUE)
med_h_founder <- median(train$h_founder, na.rm = TRUE)

fit <- readRDS(model_path)
# Optional k-predictor for recovery of demoted cases.
k_model_path <- "eval/reports/hor_model_k/ranger_model_k.rds"
fit_k <- if (file.exists(k_model_path)) readRDS(k_model_path) else NULL

d <- read.table(features_path, header = TRUE, sep = "\t",
                na.strings = "NA", stringsAsFactors = FALSE)

# Impute homology NAs.
d$h_d1[is.na(d$h_d1)] <- med_h_d1
d$h_founder[is.na(d$h_founder)] <- med_h_founder

# Predict (raw RF score).
pred <- predict(fit, data = d)
d$hor_score_raw <- pred$predictions[, "1"]

# Apply Platt calibration if available.
if (!is.null(CALIBRATION) && file.exists(CALIBRATION)) {
  cal_obj <- readRDS(CALIBRATION)
  platt_fit <- cal_obj$platt
  EPS <- 1e-4
  raw_clip <- pmin(pmax(d$hor_score_raw, EPS), 1 - EPS)
  logit_raw <- log(raw_clip / (1 - raw_clip))
  d$hor_score <- predict(platt_fit,
      newdata = data.frame(logit_score = logit_raw),
      type = "response")
} else {
  d$hor_score <- d$hor_score_raw
}

# 4-category verdict.
no_signal <- (d$s1 == 0) | (d$d1 == 0)
verdict <- character(nrow(d))
verdict[no_signal] <- "no_signal"
verdict[!no_signal & d$hor_score >= T_HIGH] <- "hor"
verdict[!no_signal & d$hor_score <  T_LOW]  <- "tandem"
verdict[!no_signal & d$hor_score >= T_LOW & d$hor_score < T_HIGH] <-
  "unresolved"

# Demote-to-unresolved refinement: when the model says HOR but the
# kite family-fit found no real founder structure (family_founder_d == 0
# or family_founder_d == family_tile_d, i.e., effectively k=1), we
# can't actually report a (founder, k, tile) hierarchy — the call is
# structurally inconsistent with HOR. Demote to unresolved.
no_family <- (d$family_founder_d == 0) |
             (d$family_tile_d == 0) |
             (d$family_founder_d == d$family_tile_d)
demoted_hor <- (verdict == "hor") & no_family

# Recovery layer (Option B): for demoted HOR cases, ask the k-predictor
# for an integer multiplicity. If d1·k_pred or d1/k_pred matches any of
# d2 / d3 within tolerance, recover as HOR with the inferred family.
n_recovered <- 0
d$k_pred <- NA_integer_
d$recovered <- FALSE
if (!is.null(fit_k) && any(demoted_hor)) {
  demoted_rows <- which(demoted_hor)
  # Predict k for demoted rows
  pred_k <- predict(fit_k, data = d[demoted_rows, ])$predictions
  k_pred_int <- pmax(2L, as.integer(round(pred_k)))
  d$k_pred[demoted_rows] <- k_pred_int

  tol_bp <- 5L
  tol_rel <- 0.02
  match_period <- function(target, candidates) {
    # Return the first candidate within ±tol_bp / ±tol_rel; else NA.
    for (c in candidates) {
      if (is.na(c) || c <= 0) next
      diff <- abs(c - target)
      if (diff <= tol_bp || diff <= tol_rel * max(c, target)) return(c)
    }
    return(NA_integer_)
  }

  for (j in seq_along(demoted_rows)) {
    i <- demoted_rows[j]
    k <- k_pred_int[j]
    d1 <- d$d1[i]
    if (d1 <= 0) next
    # Scan top-5 kite peaks (d2..d5). All are above the kite noise
    # envelope by construction. Top-5 is a balance: top-3 missed cases
    # where the founder is in d4/d5, top-10 admits weak-score peaks
    # that may be sub-period noise (TRC_1__15071287 case where d7=60
    # = 178/3 was a weak peak but matched k_pred=3).
    candidates <- c(d$d2[i], d$d3[i], d$d4[i], d$d5[i])
    # Hypothesis A: d1 is HOR-unit, candidate founder = d1/k
    cand_founder_A <- as.integer(round(d1 / k))
    matched_A <- match_period(cand_founder_A, candidates)
    # Hypothesis B: d1 is founder, candidate tile = d1*k
    cand_tile_B <- as.integer(d1 * k)
    matched_B <- match_period(cand_tile_B, candidates)

    if (!is.na(matched_A) && cand_founder_A >= 15) {
      d$family_founder_d[i] <- as.integer(matched_A)
      d$family_tile_d[i]    <- as.integer(d1)
      d$recovered[i] <- TRUE
    } else if (!is.na(matched_B) && cand_tile_B > d1) {
      d$family_founder_d[i] <- as.integer(d1)
      d$family_tile_d[i]    <- as.integer(matched_B)
      d$recovered[i] <- TRUE
    }
  }
  recovered_idx <- which(d$recovered)
  n_recovered <- length(recovered_idx)
  # Recovered cases keep verdict = "hor"; only the non-recovered demote.
  demoted_hor <- demoted_hor & !d$recovered
}
verdict[demoted_hor] <- "unresolved"
n_demoted <- sum(demoted_hor)
if (n_demoted > 0) {
  cat(sprintf("[note] demoted %d HOR→unresolved (no clean family found)\n",
              n_demoted))
}
if (n_recovered > 0) {
  cat(sprintf("[note] recovered %d HOR via k-predictor (d2/d3 matched k_pred prediction)\n",
              n_recovered))
}
d$verdict <- verdict

# For HOR: founder = family_founder_d if known, else d1; tile = family_tile_d
# if known, else NA; multiplicity = tile / founder (rounded).
d$founder <- NA_integer_
d$multiplicity <- NA_integer_
d$tile <- NA_integer_
hor_idx <- which(d$verdict == "hor")
for (i in hor_idx) {
  f <- if (d$family_founder_d[i] > 0) d$family_founder_d[i] else d$d1[i]
  t <- if (d$family_tile_d[i] > 0) d$family_tile_d[i] else d$d1[i]
  if (f > 0 && t > 0) {
    d$founder[i] <- f
    d$tile[i] <- t
    d$multiplicity[i] <- max(1L, as.integer(round(t / f)))
  }
}
# For tandem: monomer = d1, multiplicity = 1, no founder.
tan_idx <- which(d$verdict == "tandem")
d$tile[tan_idx] <- d$d1[tan_idx]
d$multiplicity[tan_idx] <- 1L

# Reorder columns for output: identifiers + verdict + score + structure +
# original features.
out_cols <- c(
  "case_id", "stratum", "array_length",
  "hor_score", "verdict", "founder", "multiplicity", "tile",
  setdiff(colnames(d), c("case_id", "stratum", "array_length",
                          "hor_score", "verdict", "founder",
                          "multiplicity", "tile"))
)
write.table(d[, out_cols], out_path, sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")

# Brief reporting summary.
cat(sprintf("wrote %d rows → %s\n", nrow(d), out_path))
cat(sprintf("thresholds: t_low=%.2f, t_high=%.2f\n", T_LOW, T_HIGH))
cat("Verdict distribution:\n")
print(table(d$verdict))
# If ground-truth labels are present (hor_signal), report metrics.
if ("hor_signal" %in% colnames(d)) {
  d$y <- as.integer(d$hor_signal >= 0.10)
  cat("\nGround-truth check (hor_signal >= 0.10):\n")
  ct <- table(verdict = d$verdict, truth = d$y)
  print(ct)
  if (sum(d$verdict == "hor") > 0) {
    tp <- sum(d$verdict == "hor" & d$y == 1)
    fp <- sum(d$verdict == "hor" & d$y == 0)
    cat(sprintf("HOR call precision: %.3f (%d/%d)\n",
                tp / (tp + fp), tp, tp + fp))
  }
  if (sum(d$verdict %in% c("tandem", "no_signal")) > 0) {
    tn <- sum(d$verdict %in% c("tandem", "no_signal") & d$y == 0)
    fn <- sum(d$verdict %in% c("tandem", "no_signal") & d$y == 1)
    cat(sprintf("NOT_HOR (tandem+no_signal) specificity: %.3f (%d/%d)\n",
                tn / (tn + fn), tn, tn + fn))
  }
}
