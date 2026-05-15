#!/usr/bin/env Rscript
# Filter ground_truth/params.tsv into the training subset:
#   - drop horsubmono / nullsubmono (structurally ambiguous — kite-valid
#     sub-motif picks). Verdict cannot be made from kite alone here.
#   - keep monomer_len ∈ [30, 3000]
#   - keep hor_order ∈ [1, 10]
#
# Output: eval/training_data/params_filtered.tsv (same schema as input).

args <- commandArgs(trailingOnly = TRUE)
in_path  <- if (length(args) >= 1) args[1] else "ground_truth/params.tsv"
out_path <- if (length(args) >= 2) args[2] else "eval/training_data/params_filtered.tsv"

p <- read.table(in_path, header = TRUE, sep = "\t",
                stringsAsFactors = FALSE, check.names = FALSE)
stratum <- sub("_[0-9]+$", "", p$case_id)
keep <- !(stratum %in% c("horsubmono", "nullsubmono")) &
        p$monomer_len >= 30 & p$monomer_len <= 3000 &
        p$hor_order >= 1 & p$hor_order <= 10

p_kept <- p[keep, , drop = FALSE]
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write.table(p_kept, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("kept %d / %d rows → %s\n", nrow(p_kept), nrow(p), out_path))
cat("Stratum distribution:\n")
print(table(sub("_[0-9]+$", "", p_kept$case_id)))
cat("hor_order distribution:\n")
print(table(p_kept$hor_order))
