#!/usr/bin/env Rscript

cat("=== Phase II – Immune Subtype Construction (TNBC-like TCGA) ===\n")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

## ---------- 1. Paths ----------

# IMPORTANT: assume working dir = /mnt/c/TNBC_project
immune_os_file <- "Immune_Spatial_Immunotherapy/Phase_I_Immune_Deconv/results/immune_scores/immune_hot_cold_TNBC_like_TCGA_OS.tsv"
tnbc_subtype_file <- "data_raw/tcga/clinical/tcga_clinical_with_subtype_proxy.csv"

out_dir <- "Immune_Spatial_Immunotherapy/Phase_II_Immune_Subtypes/results/immune_subtypes"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---------- 2. Load Phase I immune + OS table ----------

cat("Reading immune + OS table from:", immune_os_file, "\n")
imm <- fread(immune_os_file)

cat("Immune table dimensions:", paste(dim(imm), collapse = " x "), "\n")
cat("Immune columns:", paste(head(colnames(imm)), collapse = ", "), "...\n")

# Expect: submitter_id, ImmuneScore, immune_group, os_time, os_event, age_at_diagnosis, tumor_stage (some may be NA)
if (!"submitter_id" %in% names(imm)) {
  stop("submitter_id column not found in immune OS table.")
}

# Create 12-char ID
imm[, id12 := substr(submitter_id, 1, 12)]
cat("Example immune IDs (submitter_id, id12):\n")
print(head(imm[, .(submitter_id, id12)]))

## ---------- 3. Load TNBC subtype clinical file ----------

cat("Reading TNBC subtype clinical from:", tnbc_subtype_file, "\n")
clin <- fread(tnbc_subtype_file)
cat("Clinical dimensions:", paste(dim(clin), collapse = " x "), "\n")
cat("Clinical columns:", paste(colnames(clin), collapse = ", "), "\n")

# Expect columns: patient_id, ER_status, PR_status, HER2_status, TNBC, vital_status,
#                 OS_time_days, OS_status, age_at_diagnosis, ajcc_pathologic_stage, PAM50, TNBC_proxy

if (!"patient_id" %in% names(clin)) {
  stop("patient_id column not found in TNBC clinical subtype file.")
}

# Build 12-char ID
clin[, id12 := substr(patient_id, 1, 12)]
cat("Example clinical IDs (patient_id, id12):\n")
print(head(clin[, .(patient_id, id12, TNBC_proxy, PAM50)]))

## ---------- 4. Merge immune + subtype ----------

# Keep only key subtype-related columns from clinical
clin_sub <- clin[, .(
  id12,
  patient_id,
  ER_status,
  PR_status,
  HER2_status,
  TNBC,
  PAM50,
  TNBC_proxy,
  OS_time_days,
  OS_status,
  age_at_diagnosis,
  ajcc_pathologic_stage
)]

merged <- merge(
  imm,
  clin_sub,
  by = "id12",
  all.x = TRUE,
  sort = FALSE
)

cat("Merged table dimensions:", paste(dim(merged), collapse = " x "), "\n")
cat("Merged columns:", paste(colnames(merged), collapse = ", "), "\n")

na_subtype <- sum(is.na(merged$TNBC_proxy))
cat("TNBC_proxy missing in", na_subtype, "of", nrow(merged), "samples\n")

## ---------- 5. Define immune_subtype (BLIS / BLIA / IM / LAR) ----------

# Rule:
# 1) If TNBC_proxy == "LAR" → LAR
# 2) else if immune_group == "Hot" → BLIA
# 3) else if immune_group == "Cold" → BLIS
# 4) else → IM

merged[, immune_subtype := NA_character_]

merged[!is.na(TNBC_proxy) & TNBC_proxy == "LAR", immune_subtype := "LAR"]
merged[is.na(immune_subtype) & immune_group == "Hot", immune_subtype := "BLIA"]
merged[is.na(immune_subtype) & immune_group == "Cold", immune_subtype := "BLIS"]
merged[is.na(immune_subtype), immune_subtype := "IM"]

merged[, immune_subtype := factor(
  immune_subtype,
  levels = c("BLIS", "IM", "BLIA", "LAR")
)]

cat("Immune subtype counts:\n")
print(table(merged$immune_subtype, useNA = "ifany"))

## ---------- 6. Write outputs ----------

# Full merged table
full_out <- file.path(out_dir, "immune_subtypes_TNBC_like_TCGA_full.tsv")
fwrite(merged, full_out, sep = "\t")
cat("Wrote FULL immune subtype table to:", full_out, "\n")

# Core minimal table
core <- merged[, .(
  submitter_id,
  id12,
  ImmuneScore,
  immune_group,
  immune_subtype,
  TNBC_proxy,
  PAM50,
  os_time,
  os_event,
  age_at_diagnosis,
  ajcc_pathologic_stage
)]

core_out <- file.path(out_dir, "immune_subtypes_TNBC_like_TCGA_core.tsv")
fwrite(core, core_out, sep = "\t")
cat("Wrote CORE immune subtype table to:", core_out, "\n")

## ---------- 7. Barplot: subtype counts ----------

plt_counts <- ggplot(
  merged,
  aes(x = immune_subtype)
) +
  geom_bar() +
  theme_minimal(base_size = 12) +
  xlab("Immune Subtype") +
  ylab("Number of TNBC-like cases") +
  ggtitle("Immune Subtype Distribution – TNBC-like TCGA") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

png(file.path(out_dir, "immune_subtype_counts_TNBC_like_TCGA.png"),
    width = 1000, height = 800, res = 150)
print(plt_counts)
dev.off()

cat("Saved immune subtype counts plot to:",
    file.path(out_dir, "immune_subtype_counts_TNBC_like_TCGA.png"), "\n")

cat("=== DONE Phase II – Immune subtype construction ===\n")
