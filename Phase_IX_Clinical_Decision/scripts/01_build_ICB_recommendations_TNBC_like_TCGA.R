#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

cat("=== Phase IX – Clinical Immunotherapy Decision Engine (TNBC-like TCGA) ===\n\n")

# ----------------------------------------------------------------------
# 1) File paths (FIXED BASE PATH)
# ----------------------------------------------------------------------
# Working dir: /mnt/c/TNBC_project/Immune_Spatial_Immunotherapy/Phase_IX_Clinical_Decision
# So ".." = Immune_Spatial_Immunotherapy
base_dir <- normalizePath(
  file.path(getwd(), ".."),
  winslash = "/", mustWork = TRUE
)

iir_file <- file.path(
  base_dir,
  "Phase_VIII_Immunotherapy_Integration",
  "results", "integration",
  "IIR_table_TNBC_like_TCGA.tsv"
)

out_dir <- file.path(
  getwd(),
  "results", "clinical_decisions"
)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_main    <- file.path(out_dir, "ICB_recommendations_TNBC_like_TCGA.tsv")
out_summary <- file.path(out_dir, "ICB_recommendations_summary_TNBC_like_TCGA.tsv")

cat("Reading IIR table from:\n  ", iir_file, "\n")

# ----------------------------------------------------------------------
# 2) Read integrated table
# ----------------------------------------------------------------------
dt <- fread(iir_file)
cat("Loaded rows:", nrow(dt), "cols:", ncol(dt), "\n\n")

# Required columns check
required <- c(
  "submitter_id","immune_subtype","ImmuneScore","PD1_PDL1_signature",
  "Immune_Exclusion_index","Exclusion_tertile","DPI_index","DPI_tertile",
  "therapy_access_group","TMB_tertile","APOBEC_group",
  "IIR_group","IIR_score"
)

missing <- setdiff(required, colnames(dt))
if (length(missing) > 0) {
  stop("Missing columns: ", paste(missing, collapse = ", "))
}

# ----------------------------------------------------------------------
# 3) Helper biomarker flags
# ----------------------------------------------------------------------
cat("--- Deriving helper biomarkers ---\n")

# PD-L1 high/low threshold
pd1_cut <- median(dt$PD1_PDL1_signature, na.rm = TRUE)
dt[, PD1_status := ifelse(PD1_PDL1_signature >= pd1_cut, "High", "Low")]

# Immune hot flag
score_cut <- median(dt$ImmuneScore, na.rm = TRUE)

dt[, Immune_hot_flag := fifelse(
  (IIR_group == "High_ICB_ready") |
    (immune_subtype %in% c("BLIA","IM") & ImmuneScore >= score_cut),
  "Immune_hot", "Immune_cold_or_mixed"
)]

# ----------------------------------------------------------------------
# 4) Build numeric ICB_score
# ----------------------------------------------------------------------
cat("--- Computing ICB_score ---\n")

dt[, ICB_score := 0]

# IIR group
dt[IIR_group == "High_ICB_ready", ICB_score := ICB_score + 2]
dt[IIR_group == "Intermediate",   ICB_score := ICB_score + 1]

# PD1
dt[PD1_status == "High", ICB_score := ICB_score + 1]

# Subtypes
dt[immune_subtype %in% c("BLIA","IM"), ICB_score := ICB_score + 1]

# TMB & APOBEC
dt[TMB_tertile == "High",            ICB_score := ICB_score + 1]
dt[APOBEC_group == "APOBEC_high",    ICB_score := ICB_score + 1]

# Spatial penalties
dt[Exclusion_tertile == "High", ICB_score := ICB_score - 1]
dt[DPI_tertile       == "Low",  ICB_score := ICB_score - 1]

# Access
dt[therapy_access_group == "High_access", ICB_score := ICB_score + 0.5]

cat("ICB_score summary:\n")
print(summary(dt$ICB_score))
cat("\n")

# ----------------------------------------------------------------------
# 5) Convert ICB_score → priority category
# ----------------------------------------------------------------------
dt[, ICB_priority := fifelse(
  ICB_score >= 4, "High_ICB_priority",
  fifelse(ICB_score >= 2, "Intermediate_ICB_priority", "Low_ICB_priority")
)]

cat("ICB_priority counts:\n")
print(dt[, .N, by = ICB_priority])
cat("\n")

# ----------------------------------------------------------------------
# 6) Therapy recommendations (rule-based)
# ----------------------------------------------------------------------
dt[, recommended_therapy := "Undetermined"]

# High ICB priority + good spatial + good access
dt[
  ICB_priority == "High_ICB_priority" &
    therapy_access_group == "High_access" &
    Exclusion_tertile != "High",
  recommended_therapy := "ICB_monotherapy_or_ICB_plus_light_chemo"
]

# High ICB priority but bad access or high exclusion
dt[
  ICB_priority == "High_ICB_priority" &
    (therapy_access_group != "High_access" | Exclusion_tertile == "High"),
  recommended_therapy := "ICB_plus_standard_chemo"
]

# Intermediate
dt[
  ICB_priority == "Intermediate_ICB_priority",
  recommended_therapy := "ICB_plus_chemo_or_combo_strategies"
]

# Low
dt[
  ICB_priority == "Low_ICB_priority",
  recommended_therapy := "Non_ICB_chemo_or_targeted_plus_TME_modulation"
]

# ----------------------------------------------------------------------
# 7) Rationale string
# ----------------------------------------------------------------------
dt[, rationale_short := paste0(
  "Subtype=", immune_subtype,
  "; IIR=", IIR_group,
  "; PD1=", PD1_status,
  "; TMB=", TMB_tertile,
  "; APOBEC=", APOBEC_group,
  "; Exclusion=", Exclusion_tertile,
  "; DPI=", DPI_tertile,
  "; Access=", therapy_access_group
)]

# ----------------------------------------------------------------------
# 8) Save main output
# ----------------------------------------------------------------------
cat("--- Writing output tables ---\n")

fwrite(dt[, .(
  submitter_id, IIR_group, IIR_score,
  ImmuneScore, immune_subtype, PD1_status,
  TMB_tertile, APOBEC_group,
  Exclusion_tertile, DPI_tertile,
  therapy_access_group,
  ICB_score, ICB_priority,
  recommended_therapy, rationale_short
)], out_main, sep = "\t")

fwrite(dt[, .N, by = recommended_therapy], out_summary, sep = "\t")

cat("Saved:\n  ", out_main, "\n  ", out_summary, "\n\n")
cat("=== DONE Phase IX – Decision Engine ===\n")