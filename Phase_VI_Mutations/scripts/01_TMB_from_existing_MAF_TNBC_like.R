#!/usr/bin/env Rscript

cat("=== Phase VI – TMB from existing MAF (TNBC-like TCGA) ===\n")

suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------------------------------------------
# 1) Locate an existing MAF file in your TNBC project
# ------------------------------------------------------------------
project_root <- normalizePath(file.path("..", ".."), mustWork = FALSE)
cat("Searching for .maf files under:\n  ", project_root, "\n")

maf_candidates <- list.files(
  path       = project_root,
  pattern    = "\\.maf$",
  full.names = TRUE,
  recursive  = TRUE
)

if (length(maf_candidates) == 0) {
  stop(
    "ERROR: No .maf files found under ", project_root, "\n",
    "If you know the MAF path (e.g. data_proc/TNBC_only.maf), tell me and we will hard-code it."
  )
}

cat("MAF candidates found:\n")
print(maf_candidates)

# Prefer TNBC-only / BRCA / mutect files if present
preferred <- maf_candidates[grepl("TNBC_only|BRCA|mutect", maf_candidates, ignore.case = TRUE)]
if (length(preferred) > 0) {
  maf_file <- preferred[1]
} else {
  maf_file <- maf_candidates[1]
}

cat("Using MAF file:\n  ", maf_file, "\n")

# ------------------------------------------------------------------
# 2) Load MAF and extract non-silent mutations
# ------------------------------------------------------------------
maf <- fread(maf_file)

col_barcode <- intersect(
  c("Tumor_Sample_Barcode", "Tumor.Sample.Barcode", "Sample"),
  names(maf)
)
col_class <- intersect(
  c("Variant_Classification", "Variant.Classification"),
  names(maf)
)

if (length(col_barcode) == 0 || length(col_class) == 0) {
  stop(
    "ERROR: Could not find Tumor_Sample_Barcode / Variant_Classification columns in MAF.\n",
    "Found columns:\n", paste(names(maf), collapse = ", ")
  )
}

setnames(maf, col_barcode[1], "Tumor_Sample_Barcode")
setnames(maf, col_class[1],   "Variant_Classification")

non_silent <- c(
  "Frame_Shift_Del", "Frame_Shift_Ins",
  "Missense_Mutation", "Nonsense_Mutation",
  "Nonstop_Mutation", "Splice_Site",
  "Translation_Start_Site",
  "In_Frame_Del", "In_Frame_Ins"
)

maf_ns <- maf[Variant_Classification %in% non_silent]

cat("Total variants in MAF:      ", nrow(maf), "\n")
cat("Non-silent variants retained:", nrow(maf_ns), "\n")

# Case-level (12-char) IDs
maf_ns[, submitter_id := substr(Tumor_Sample_Barcode, 1, 12)]

tmb_counts <- maf_ns[, .(n_nonsilent = .N), by = submitter_id]

# Approximate exome size in Mb -> just for scaling, relative TMB is what matters
exome_mb <- 38
tmb_counts[, TMB := n_nonsilent / exome_mb]

cat("Example TMB rows:\n")
print(head(tmb_counts))

# ------------------------------------------------------------------
# 3) Merge with TNBC-like immune OS table
# ------------------------------------------------------------------
immune_file <- "../Phase_I_Immune_Deconv/results/immune_scores/immune_hot_cold_TNBC_like_TCGA_OS.tsv"
cat("Reading immune OS file from:\n  ", immune_file, "\n")

immune <- fread(immune_file)

merged <- merge(
  immune,
  tmb_counts,
  by      = "submitter_id",
  all.x   = TRUE
)

cat("Rows in immune OS:   ", nrow(immune), "\n")
cat("Rows with TMB values:", sum(!is.na(merged$TMB)), "\n")

# ------------------------------------------------------------------
# 4) Save outputs
#    (a) Phase VI mutations folder (rich table)
#    (b) Phase III immunotherapy folder (minimal TMB table)
# ------------------------------------------------------------------
out_dir_phase <- "results/mutations"
if (!dir.exists(out_dir_phase)) {
  dir.create(out_dir_phase, recursive = TRUE, showWarnings = FALSE)
}

out_file_phase <- file.path(out_dir_phase, "TMB_TNBC_like_TCGA_from_MAF.tsv")

# IMPORTANT: only use columns that actually exist in immune_OS + TMB
keep_cols_phase <- intersect(
  c("submitter_id",
    "ImmuneScore",
    "immune_group",
    "os_time",
    "os_event",
    "age_at_diagnosis",
    "tumor_stage",
    "n_nonsilent",
    "TMB"),
  names(merged)
)

fwrite(
  merged[, ..keep_cols_phase],
  out_file_phase,
  sep = "\t"
)

cat("Saved Phase VI TMB table to:\n  ", out_file_phase, "\n")

# Also export minimal file for Phase III ICB composite models
out_dir_icb <- "../Phase_III_Immunotherapy_Prediction/results/immunotherapy"
if (!dir.exists(out_dir_icb)) {
  dir.create(out_dir_icb, recursive = TRUE, showWarnings = FALSE)
}

out_file_icb <- file.path(out_dir_icb, "TMB_TNBC_like_TCGA.tsv")

fwrite(
  merged[, .(submitter_id, n_nonsilent, TMB)],
  out_file_icb,
  sep = "\t"
)

cat("Saved TMB file for ICB models to:\n  ", out_file_icb, "\n")
cat("=== DONE – TMB from existing MAF ===\n")
