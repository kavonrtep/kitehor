#!/usr/bin/env Rscript
# Export ranger random-forest models to a compact JSON format that the
# Rust classifier can load and traverse.
#
# Output schema (per model file):
#   {
#     "treetype": "probability"|"regression",
#     "num_trees": int,
#     "feature_names": [str, ...],
#     "class_levels": [str, ...]                  (probability mode only)
#     "trees": [
#       {
#         "split_var":   [int, ...]    (0-indexed; -1 = leaf)
#         "split_val":   [f64, ...]    (0.0 at leaves)
#         "left":        [u32, ...]    (0 at leaves)
#         "right":       [u32, ...]    (0 at leaves)
#         "prediction":  [f64, ...]    (P(class=level2) for prob,
#                                       mean for regression; 0 at internals)
#       }, ...
#     ]
#   }
#
# Split semantics (ranger): for numeric features, x <= split_val goes
# to "left" (verified empirically). All features in our models are
# numeric, so categorical splits are not handled here.
#
# Usage:
#   Rscript export_ranger.R --in <model.rds> --out <model.json>

suppressPackageStartupMessages({ library(ranger) })

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(key, default = NULL) {
  i <- which(args == key)
  if (length(i) == 0) return(default)
  args[i + 1]
}
in_path  <- parse_arg("--in")
out_path <- parse_arg("--out")
if (is.null(in_path) || is.null(out_path)) {
  stop("usage: export_ranger.R --in <model.rds> --out <model.json>")
}

fit <- readRDS(in_path)
ttype <- fit$treetype
if (!(ttype %in% c("Probability estimation", "Regression"))) {
  stop(sprintf("unsupported treetype: %s", ttype))
}
is_prob <- ttype == "Probability estimation"
cat(sprintf("input: %s\n", in_path))
cat(sprintf("treetype: %s, num.trees: %d, num.features: %d\n",
            ttype, fit$num.trees, fit$num.independent.variables))

feat_names <- fit$forest$independent.variable.names
class_levels <- NULL
if (is_prob) {
  class_levels <- fit$forest$levels
  if (is.null(class_levels)) {
    class_levels <- as.character(fit$forest$class.values)
  }
}

# Hand-rolled JSON writer. We avoid jsonlite (not always available) and
# emit a strict-but-compact JSON document.
quote_str <- function(s) {
  s <- gsub("\\\\", "\\\\\\\\", s)
  s <- gsub("\"", "\\\\\"", s)
  paste0("\"", s, "\"")
}
fmt_int_arr <- function(v) paste0("[", paste(as.integer(v), collapse=","), "]")
fmt_num_arr <- function(v) {
  # Use %.17g to round-trip every IEEE-754 double exactly.
  paste0("[", paste(sprintf("%.17g", v), collapse=","), "]")
}
fmt_str_arr <- function(v) paste0("[", paste(vapply(v, quote_str, character(1)), collapse=","), "]")

fh <- file(out_path, open = "w")
on.exit(close(fh), add = TRUE)
writeLines("{", fh)
writeLines(sprintf("  \"treetype\": %s,",
                   quote_str(if (is_prob) "probability" else "regression")), fh)
writeLines(sprintf("  \"num_trees\": %d,", fit$num.trees), fh)
writeLines(sprintf("  \"feature_names\": %s,", fmt_str_arr(feat_names)), fh)
if (is_prob) {
  writeLines(sprintf("  \"class_levels\": %s,", fmt_str_arr(class_levels)), fh)
}
writeLines("  \"trees\": [", fh)

for (t in seq_len(fit$num.trees)) {
  ti <- treeInfo(fit, tree = t)
  n <- nrow(ti)
  split_var  <- integer(n)
  split_val  <- numeric(n)
  left_arr   <- integer(n)
  right_arr  <- integer(n)
  pred       <- numeric(n)
  term <- ti$terminal
  internal <- !term
  split_var[term]      <- -1L
  split_var[internal]  <- as.integer(ti$splitvarID[internal])
  split_val[internal]  <- as.numeric(ti$splitval[internal])
  left_arr[internal]   <- as.integer(ti$leftChild[internal])
  right_arr[internal]  <- as.integer(ti$rightChild[internal])
  if (is_prob) {
    pred[term] <- ti$pred.1[term]   # P(class = second level) = HOR
  } else {
    pred[term] <- ti$prediction[term]
  }
  sep <- if (t < fit$num.trees) "," else ""
  writeLines(sprintf("    {\"split_var\":%s,\"split_val\":%s,\"left\":%s,\"right\":%s,\"prediction\":%s}%s",
                     fmt_int_arr(split_var),
                     fmt_num_arr(split_val),
                     fmt_int_arr(left_arr),
                     fmt_int_arr(right_arr),
                     fmt_num_arr(pred),
                     sep), fh)
}
writeLines("  ]", fh)
writeLines("}", fh)

cat(sprintf("wrote: %s (%d trees)\n", out_path, fit$num.trees))
