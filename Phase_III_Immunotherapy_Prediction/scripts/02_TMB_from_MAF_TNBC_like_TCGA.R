#!/usr/bin/env Rscript

cat("=== Phase III – TMB from MAF (TNBC-like TCGA) ===\n")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

maf_file       <- "inputs/mutations/tcga_brca_mc3_public.maf"
immune_os_file <- "../Phase_I_Immune_Deconv/results/immune_scores/immune_hot_cold_TNBC_like_TCGA_OS.tsv"
out_dir        <- "results/immunotherapy"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Using immune OS file:", immune_os_file, "\n")
immune_os <- fread(immune_os_file)

if (!"submitter_id" %in% names(immune_os)) {
  stop("submitter_id column not found in immune OS file: ", immune_os_file)
}

## ------------------------------------------------------------------
## CASE 1: MAF file is missing  → create NA TMB table (graceful path)
## ------------------------------------------------------------------
if (!file.exists(maf_file)) {
  cat("WARNING: MAF file not found at:\n  ", maf_file, "\n", sep = "")
  cat("→ Creating placeholder TMB table with NA values so downstream scripts still work.\n")

  tmb_tnbc <- immune_os %>%
    select(submitter_id) %>%
    mutate(
      n_nonsilent = NA_integer_,
      TMB         = NA_real_
    )

  out_file <- file.path(out_dir, "TMB_TNBC_like_TCGA.tsv")
  fwrite(as.data.table(tmb_tnbc), out_file, sep = "\t")
  cat("Wrote placeholder TMB table with NA values to:\n  ", out_file, "\n", sep = "")
  cat("=== DONE TMB Phase (placeholder, no MAF) ===\n")
  quit(save = "no")
}

## ------------------------------------------------------------------
## CASE 2: MAF file exists  → compute real TMB
## ------------------------------------------------------------------
cat("Reading MAF from:", maf_file, "\n")
maf <- fread(maf_file)

required_cols <- c("Tumor_Sample_Barcode", "Variant_Classification")
missing_cols  <- setdiff(required_cols, names(maf))
if (length(missing_cols) > 0) {
  stop("MAF missing required columns: ", paste(missing_cols, collapse = ", "))
}

cat("Total MAF rows:", nrow(maf), "\n")

nonsilent_classes <- c(
  "Frame_Shift_Del", "Frame_Shift_Ins", "Missense_Mutation",
  "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site",
  "Translation_Start_Site", "In_Frame_Ins", "In_Frame_Del"
)

maf_ns <- maf[Variant_Classification %in% nonsilent_classes]
cat("Non-silent variants:", nrow(maf_ns), "\n")

maf_ns[, case_id := substr(Tumor_Sample_Barcode, 1, 12)]
tmb_counts <- maf_ns[, .N, by = case_id]
setnames(tmb_counts, "N", "n_nonsilent")

exome_mb <- 38
tmb_counts[, TMB := n_nonsilent / exome_mb]

cat("TMB table rows:", nrow(tmb_counts), "\n")
cat("TMB summary:\n")
print(summary(tmb_counts$TMB))

tmb_tnbc <- immune_os %>%
  select(submitter_id) %>%
  left_join(
    as.data.frame(tmb_counts),
    by = c("submitter_id" = "case_id")
  )

cat("Rows after restricting to TNBC-like:", nrow(tmb_tnbc), "\n")
cat("Non-NA TMB rows:", sum(!is.na(tmb_tnbc$TMB)), "\n")

out_file <- file.path(out_dir, "TMB_TNBC_like_TCGA.tsv")
fwrite(as.data.table(tmb_tnbc), out_file, sep = "\t")
cat("Wrote TMB table to:", out_file, "\n")

cat("=== DONE TMB Phase (real MAF) ===\n")
