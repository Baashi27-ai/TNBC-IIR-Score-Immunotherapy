#!/usr/bin/env Rscript

cat("=== Phase XI – Escape-Adjusted ICB Score (TNBC-like TCGA) ===\n")

suppressPackageStartupMessages({
  library(data.table)
})

## 1) INPUT: HLA escape–merged table from Phase X
hla_file <- "../Phase_X_HLA_ImmuneEscape/results/hla/HLA_escape_TNBC_like_TCGA.tsv"

if (!file.exists(hla_file)) {
  stop("HLA escape file not found at: ", hla_file)
}

dt <- fread(hla_file)
cat("Loaded HLA escape table:\n  ", hla_file, "\n")
cat("Dimensions:", nrow(dt), "rows x", ncol(dt), "cols\n")

cat("Columns:\n")
print(colnames(dt))

## 2) Basic sanity checks
required_cols <- c(
  "submitter_id", "os_time", "os_event",
  "ImmuneScore", "immune_subtype",
  "IIR_group", "IIR_score",
  "PD1_PDL1_signature",
  "TMB", "TMB_tertile",
  "APOBEC_group",
  "HLA_escape_group", "MHC_I_axis_tertile",
  "therapy_access_group"
)

missing <- setdiff(required_cols, colnames(dt))
if (length(missing) > 0) {
  stop("Missing required columns in HLA escape table: ",
       paste(missing, collapse = ", "))
}

## 3) Derive PD1 status (High / Low) based on median
dt[, PD1_PDL1_signature := as.numeric(PD1_PDL1_signature)]
pd1_med <- median(dt$PD1_PDL1_signature, na.rm = TRUE)
cat("PD1/PD-L1 signature median:", pd1_med, "\n")

dt[, PD1_status := fifelse(
  is.na(PD1_PDL1_signature),
  NA_character_,
  fifelse(PD1_PDL1_signature >= pd1_med, "High", "Low")
)]

## 4) Scoring components

# Base from IIR_group
# High_ICB_ready = 2, Intermediate = 1, Poor_ICB_ready = 0
dt[, IIR_base_score := fifelse(
  IIR_group == "High_ICB_ready", 2,
  fifelse(IIR_group == "Intermediate", 1,
          fifelse(IIR_group == "Poor_ICB_ready", 0, NA_real_))
)]

# PD1 contribution: High = +0.5, Low = 0
dt[, PD1_score := fifelse(PD1_status == "High", 0.5,
                          fifelse(PD1_status == "Low", 0, NA_real_))]

# TMB contribution from tertiles: High=+1, Mid=+0.5, Low=0
dt[, TMB_score := fifelse(
  TMB_tertile == "High", 1,
  fifelse(TMB_tertile == "Mid", 0.5,
          fifelse(TMB_tertile == "Low", 0, NA_real_))
)]

# APOBEC contribution: APOBEC_high=+0.5, APOBEC_low=0
dt[, APOBEC_score := fifelse(
  APOBEC_group == "APOBEC_high", 0.5,
  fifelse(APOBEC_group == "APOBEC_low", 0, NA_real_)
)]

# Therapy access contribution: High_access=+0.5, Medium=0, Poor_access=-0.5
dt[, Access_score := fifelse(
  therapy_access_group == "High_access", 0.5,
  fifelse(therapy_access_group == "Medium", 0,
          fifelse(therapy_access_group == "Poor_access", -0.5, NA_real_))
)]

# HLA immune escape penalty:
# No_escape = 0, Partial_escape = -1, (if ever Strong_escape) = -2
dt[, Escape_penalty := fifelse(
  HLA_escape_group == "No_escape", 0,
  fifelse(HLA_escape_group == "Partial_escape", -1,
          fifelse(HLA_escape_group == "Strong_escape", -2, NA_real_))
)]

## 5) Final escape-adjusted ICB score
dt[, escape_adjusted_ICB_score :=
      IIR_base_score +
      PD1_score +
      TMB_score +
      APOBEC_score +
      Access_score +
      Escape_penalty]

cat("Escape-adjusted ICB score summary:\n")
print(summary(dt$escape_adjusted_ICB_score))

## 6) Group classification
# Rough rule-of-thumb:
#   >= 3.5       → Escape_adjusted_High_ICB_candidate
#   2.0–<3.5     → Escape_adjusted_Intermediate_candidate
#   < 2.0        → Escape_adjusted_Poor_ICB_candidate

dt[, escape_adjusted_ICB_group := fifelse(
  is.na(escape_adjusted_ICB_score), NA_character_,
  fifelse(escape_adjusted_ICB_score >= 3.5, "Escape_adjusted_High_ICB_candidate",
          fifelse(escape_adjusted_ICB_score >= 2.0, "Escape_adjusted_Intermediate_candidate",
                  "Escape_adjusted_Poor_ICB_candidate"))
)]

cat("Escape-adjusted ICB group counts:\n")
print(dt[, .N, by = escape_adjusted_ICB_group][order(escape_adjusted_ICB_group)])

## 7) Output minimal + key columns
out_dir  <- "results/escape_adjusted"
out_file <- file.path(out_dir, "escape_adjusted_ICB_TNBC_like_TCGA.tsv")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

keep_cols <- c(
  "submitter_id", "id12",
  "ImmuneScore", "immune_group", "immune_subtype",
  "TNBC_proxy", "PAM50",
  "os_time", "os_event",
  "age_at_diagnosis", "ajcc_pathologic_stage",
  "PD1_PDL1_signature", "PD1_status",
  "TMB", "TMB_tertile",
  "APOBEC_group",
  "Immune_composite", "Stroma_composite",
  "Immune_Stroma_ratio", "Immune_Exclusion_index",
  "Immune_Stroma_tertile",
  "DPI_index", "DPI_tertile",
  "Exclusion_tertile", "therapy_access_group",
  "n_nonsilent", "total_snvs",
  "Ageing_CpG_prop", "APOBEC_prop",
  "IIR_score", "IIR_group", "IIR_tertile",
  "HLA_LOH_proxy_group", "MHC_I_axis", "MHC_I_axis_tertile",
  "HLA_escape_score", "HLA_escape_group",
  "IIR_base_score", "PD1_score", "TMB_score",
  "APOBEC_score", "Access_score", "Escape_penalty",
  "escape_adjusted_ICB_score", "escape_adjusted_ICB_group"
)

keep_cols <- intersect(keep_cols, colnames(dt))
dt_out <- dt[, ..keep_cols]

fwrite(dt_out, out_file, sep = "\t")
cat("Saved escape-adjusted ICB table to:\n  ", out_file, "\n")

cat("Preview:\n")
print(head(dt_out, 6))

cat("=== DONE Phase XI – Escape-Adjusted ICB Score ===\n")