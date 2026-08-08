#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

cat("=== Phase X – Merge LOH + MHC-I Axis + IIR (TNBC-like TCGA) ===\n")

setwd("/mnt/c/TNBC_project")

## 1) File paths
loh_file <- "Immune_Spatial_Immunotherapy/Phase_X_HLA_LOH/results/hla_loh/HLA_LOH_proxy_TNBC_like_TCGA.tsv"
mhc_file <- "Immune_Spatial_Immunotherapy/Phase_X_HLA_ImmuneEscape/results/hla/HLA_expression_axis_TNBC_like_TCGA.tsv"
iir_file <- "Immune_Spatial_Immunotherapy/Phase_VIII_Immunotherapy_Integration/results/integration/IIR_table_TNBC_like_TCGA.tsv"

out_dir  <- "Immune_Spatial_Immunotherapy/Phase_X_HLA_ImmuneEscape/results/hla"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "HLA_escape_TNBC_like_TCGA.tsv")

## 2) Read inputs
cat("Reading LOH proxy from:\n  ", loh_file, "\n")
loh <- fread(loh_file)

cat("Reading MHC-I expression axis from:\n  ", mhc_file, "\n")
mhc <- fread(mhc_file)

cat("Reading IIR table from:\n  ", iir_file, "\n")
iir <- fread(iir_file)

if (!"submitter_id" %in% names(loh)) stop("LOH table missing 'submitter_id'.")
if (!"submitter_id" %in% names(mhc)) stop("MHC table missing 'submitter_id'.")
if (!"submitter_id" %in% names(iir)) stop("IIR table missing 'submitter_id'.")

cat("LOH rows:", nrow(loh), "\n")
cat("MHC rows:", nrow(mhc), "\n")
cat("IIR rows:", nrow(iir), "\n")

## 3) Reduce LOH + MHC to relevant columns
loh_keep_cols <- intersect(
  c("submitter_id", "HLA_total_variants", "HLA_damaging_variants",
    "HLA_damaging_prop", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
    "HLA_LOH_proxy_group"),
  names(loh)
)

loh_small <- loh[, ..loh_keep_cols]

mhc_keep_cols <- intersect(
  c("submitter_id", "MHC_I_axis", "MHC_I_axis_tertile"),
  names(mhc)
)
mhc_small <- mhc[, ..mhc_keep_cols]

## 4) Merge: IIR is main backbone (215 rows)
dt <- merge(iir, loh_small, by = "submitter_id", all.x = TRUE)
dt <- merge(dt, mhc_small, by = "submitter_id", all.x = TRUE)

cat("Merged rows (IIR backbone):", nrow(dt), "\n")

## 5) Define HLA escape flags and score
dt[
  ,
  HLA_LOH_flag := fifelse(
    !is.na(HLA_LOH_proxy_group) & HLA_LOH_proxy_group == "High",
    1L, 0L
  )
]

dt[
  ,
  MHC_I_low_flag := fifelse(
    !is.na(MHC_I_axis_tertile) & MHC_I_axis_tertile == "Low",
    1L, 0L
  )
]

dt[, HLA_escape_score := HLA_LOH_flag + MHC_I_low_flag]

dt[
  ,
  HLA_escape_group := fifelse(
    is.na(HLA_escape_score), NA_character_,
    fifelse(
      HLA_escape_score == 0, "No_escape",
      fifelse(
        HLA_escape_score == 1, "Partial_escape",
        "Strong_escape"
      )
    )
  )
]

## 6) Save and quick summary
fwrite(dt, out_file, sep = "\t")

cat("Saved HLA escape table to:\n  ", out_file, "\n")

cat("HLA_escape_group counts:\n")
print(dt[, .N, by = HLA_escape_group])

cat("Preview:\n")
print(head(dt[
  ,
  .(submitter_id,
    IIR_group,
    IIR_score,
    HLA_LOH_proxy_group,
    MHC_I_axis_tertile,
    HLA_escape_score,
    HLA_escape_group)
]))

cat("=== DONE Phase X – HLA Immune Escape Integration ===\n")