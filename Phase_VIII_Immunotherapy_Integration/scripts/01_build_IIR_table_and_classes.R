#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

cat("=== Phase VIII – Build Integrated Immunotherapy Readiness (IIR) Table ===\n")

## ----------------------------
## 1) Input paths
## ----------------------------

therapy_file <- "../Phase_V_Drug_Penetration/results/drug_penetration/therapy_accessibility_TNBC_like_TCGA.tsv"
tmb_file     <- "../Phase_III_Immunotherapy_Prediction/results/immunotherapy/TMB_TNBC_like_TCGA.tsv"

# You confirmed this file exists:
#   C:\\TNBC_project\\Immune_Spatial_Immunotherapy\\Phase_VII_Mutational_Signatures\\results\\signatures\\signature_scores_TNBC_like_TCGA.tsv
# WSL path:
sig_file     <- "../Phase_VII_Mutational_Signatures/results/signatures/signature_scores_TNBC_like_TCGA.tsv"

out_dir  <- "results/integration"
out_file <- file.path(out_dir, "IIR_table_TNBC_like_TCGA.tsv")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

## ----------------------------
## 2) Load therapy accessibility table
## ----------------------------

cat("Reading base therapy accessibility table from:\n  ", therapy_file, "\n")
therapy <- fread(therapy_file)

cat("Base table dimensions:", nrow(therapy), "rows x", ncol(therapy), "cols\n")
cat("Base columns:\n")
print(colnames(therapy))

## ----------------------------
## 3) Load TMB table
## ----------------------------

cat("Reading TMB table from:\n  ", tmb_file, "\n")
tmb <- fread(tmb_file)

cat("TMB table columns:\n")
print(colnames(tmb))

# Keep only TMB-related columns for merge
tmb_sub <- tmb[, .(submitter_id, n_nonsilent, TMB)]

## ----------------------------
## 4) Load mutational signature summary
##     (we use signature_scores_TNBC_like_TCGA.tsv directly)
## ----------------------------

cat("Reading mutational signature summary (scores) from:\n  ", sig_file, "\n")
sig_dt <- fread(sig_file)

cat("Signature table dimensions:", nrow(sig_dt), "rows x", ncol(sig_dt), "cols\n")
cat("Signature columns:\n")
print(colnames(sig_dt))
# Expected: submitter_id, total_snvs, Ageing_CpG_prop, APOBEC_prop

## ----------------------------
## 5) Merge all layers
## ----------------------------

# Merge therapy + TMB
dt <- merge(
  therapy,
  tmb_sub,
  by = "submitter_id",
  all.x = TRUE,
  suffixes = c("", "_TMB")
)

cat("After merging TMB: ", nrow(dt), "rows x", ncol(dt), "cols\n")

# Merge signatures
dt <- merge(
  dt,
  sig_dt,
  by = "submitter_id",
  all.x = TRUE,
  suffixes = c("", "_SIG")
)

cat("After merging signatures: ", nrow(dt), "rows x", ncol(dt), "cols\n")

## ----------------------------
## 6) Quick OS / event sanity check
## ----------------------------

if (all(c("os_time", "os_event") %in% colnames(dt))) {
  cat("OS summary:\n")
  print(summary(dt$os_time))
  cat("Events (os_event == 1):", sum(dt$os_event == 1, na.rm = TRUE), "\n")
} else {
  cat("WARNING: os_time / os_event not found in merged table.\n")
}

## ----------------------------
## 7) Helper: z-score scaling
## ----------------------------

zscale <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  m <- mean(x, na.rm = TRUE)
  s <- sd(x,   na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - m) / s
}

## ----------------------------
## 8) Create TMB / APOBEC categories
## ----------------------------

# TMB tertiles (within TNBC-like cohort)
if ("TMB" %in% colnames(dt)) {
  q <- quantile(dt$TMB, probs = c(1/3, 2/3), na.rm = TRUE)
  dt[, TMB_tertile := fifelse(
    is.na(TMB), NA_character_,
    fifelse(TMB <= q[1], "Low",
      fifelse(TMB >= q[2], "High", "Mid")
    )
  )]
} else {
  dt[, TMB_tertile := NA_character_]
}

# APOBEC high vs low (median split)
if ("APOBEC_prop" %in% colnames(dt)) {
  med_apobec <- median(dt$APOBEC_prop, na.rm = TRUE)
  dt[, APOBEC_group := fifelse(
    is.na(APOBEC_prop), NA_character_,
    fifelse(APOBEC_prop >= med_apobec, "APOBEC_high", "APOBEC_low")
  )]
} else {
  dt[, APOBEC_group := NA_character_]
}

## ----------------------------
## 9) Construct Integrated Immunotherapy Readiness (IIR) score
##     Components:
##       + ImmuneScore (higher = more immune inflamed)
##       + PD1_PDL1_signature (checkpoint activation)
##       + TMB (more mutations → more neoantigens)
##       + APOBEC_prop (mutational activity often linked to ICB response)
##       - Immune_Exclusion_index (exclusion = bad)
##       + DPI_index (drug penetration)
## ----------------------------

needed_cols <- c(
  "ImmuneScore",
  "PD1_PDL1_signature",
  "TMB",
  "APOBEC_prop",
  "Immune_Exclusion_index",
  "DPI_index"
)

cat("Checking availability of IIR components:\n")
print(needed_cols)
print(colnames(dt)[colnames(dt) %in% needed_cols])

# Build scaled components (only for those that exist)
if (!("ImmuneScore" %in% colnames(dt))) dt[, ImmuneScore := NA_real_]
if (!("PD1_PDL1_signature" %in% colnames(dt))) dt[, PD1_PDL1_signature := NA_real_]
if (!("TMB" %in% colnames(dt))) dt[, TMB := NA_real_]
if (!("APOBEC_prop" %in% colnames(dt))) dt[, APOBEC_prop := NA_real_]
if (!("Immune_Exclusion_index" %in% colnames(dt))) dt[, Immune_Exclusion_index := NA_real_]
if (!("DPI_index" %in% colnames(dt))) dt[, DPI_index := NA_real_]

dt[, IIR_ImmuneScore_z      := zscale(ImmuneScore)]
dt[, IIR_PD1sig_z           := zscale(PD1_PDL1_signature)]
dt[, IIR_TMB_z              := zscale(TMB)]
dt[, IIR_APOBEC_z           := zscale(APOBEC_prop)]
dt[, IIR_Exclusion_minus_z  := zscale(-Immune_Exclusion_index)]  # negative because higher exclusion is bad
dt[, IIR_DPI_z              := zscale(DPI_index)]

# Row-wise average of available components
dt[, IIR_score := rowMeans(
  .SD,
  na.rm = TRUE
), .SDcols = c(
  "IIR_ImmuneScore_z",
  "IIR_PD1sig_z",
  "IIR_TMB_z",
  "IIR_APOBEC_z",
  "IIR_Exclusion_minus_z",
  "IIR_DPI_z"
)]

## ----------------------------
## 10) IIR tertiles & groups
## ----------------------------

valid_iir <- !is.na(dt$IIR_score)
if (any(valid_iir)) {
  q_iir <- quantile(dt$IIR_score[valid_iir], probs = c(1/3, 2/3), na.rm = TRUE)

  dt[valid_iir, IIR_tertile := fifelse(
    IIR_score <= q_iir[1], "Low",
    fifelse(IIR_score >= q_iir[2], "High", "Mid")
  )]

  dt[!valid_iir, IIR_tertile := NA_character_]

  dt[, IIR_group := factor(
    fifelse(
      is.na(IIR_tertile), NA_character_,
      fifelse(
        IIR_tertile == "High", "High_ICB_ready",
        fifelse(IIR_tertile == "Mid", "Intermediate", "Poor_ICB_ready")
      )
    ),
    levels = c("Poor_ICB_ready", "Intermediate", "High_ICB_ready")
  )]
} else {
  dt[, IIR_tertile := NA_character_]
  dt[, IIR_group   := factor(NA_character_)]
}

cat("IIR group counts:\n")
print(table(dt$IIR_group, useNA = "ifany"))

## ----------------------------
## 11) Save output
## ----------------------------

cat("Writing IIR table to:\n  ", out_file, "\n")
fwrite(dt, out_file, sep = "\t")

cat("=== DONE building IIR table ===\n")