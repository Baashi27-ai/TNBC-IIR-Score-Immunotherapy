#!/usr/bin/env Rscript

message("=== Phase I – Immune Deconvolution Summary (TNBC-like TCGA) ===")

suppressPackageStartupMessages({
  library(data.table)
})

xcell_path <- file.path("results", "immune_scores", "xcell_scores_TNBC_like_TCGA.csv")
clin_path  <- file.path("inputs", "clinical", "tcga_brca_clinical_case12.csv")
out_dir    <- file.path("results", "immune_scores")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---------- 1) Read xCell scores ----------
message("Reading xCell scores from: ", xcell_path)
xc_raw <- fread(xcell_path)
# First column is cell type name, others are samples (case IDs)
cell_types <- xc_raw[[1]]
sample_ids <- colnames(xc_raw)[-1]

xc_mat <- as.matrix(xc_raw[, -1])
rownames(xc_mat) <- cell_types
colnames(xc_mat) <- sample_ids

message("xCell matrix: ", nrow(xc_mat), " cell types x ", ncol(xc_mat), " samples")

# Transpose to samples x cell types
xc_t <- t(xc_mat)
# Convert to data.table
xc_dt <- as.data.table(xc_t, keep.rownames = "submitter_id")

# Save patient-wise table
out_samples <- file.path(out_dir, "xcell_scores_TNBC_like_TCGA_samples.tsv")
message("Writing sample-wise xCell table to: ", out_samples)
fwrite(xc_dt, out_samples, sep = "\t")

# ---------- 2) Define an ImmuneScore & hot/cold groups ----------
# xCell usually has an 'ImmuneScore' row. If missing, use mean of all rows as fallback.
if ("ImmuneScore" %in% rownames(xc_mat)) {
  immune_score <- as.numeric(xc_mat["ImmuneScore", ])
  names(immune_score) <- colnames(xc_mat)
  message("Using xCell ImmuneScore row.")
} else {
  immune_score <- colMeans(xc_mat)
  names(immune_score) <- colnames(xc_mat)
  message("No 'ImmuneScore' row found; using mean of all xCell cell scores per sample.")
}

immune_dt <- data.table(
  submitter_id = names(immune_score),
  ImmuneScore  = as.numeric(immune_score)
)

# Tertiles: Cold (bottom 1/3), Intermediate (middle), Hot (top 1/3)
q <- quantile(immune_dt$ImmuneScore, probs = c(1/3, 2/3), na.rm = TRUE)

immune_dt[, immune_group := fifelse(
  ImmuneScore <= q[1], "Cold",
  fifelse(ImmuneScore >= q[2], "Hot", "Intermediate")
)]

message("ImmuneScore summary:")
print(summary(immune_dt$ImmuneScore))
message("Group counts:")
print(table(immune_dt$immune_group))

# ---------- 3) Merge with TCGA clinical for OS ----------
message("Reading TCGA clinical file from: ", clin_path)
clin <- fread(clin_path)

# Standardize ID column used earlier ('submitter_id')
if (!"submitter_id" %in% colnames(clin)) {
  stop("Clinical file has no 'submitter_id' column – please check tcga_brca_clinical_case12.csv.")
}

clin[, submitter_id := gsub("\\.", "-", submitter_id)]

# Create OS time & event from days_to_death / days_to_last_follow_up + vital_status
if (!all(c("vital_status", "days_to_death", "days_to_last_follow_up") %in% colnames(clin))) {
  stop("Clinical file is missing vital_status/days_to_death/days_to_last_follow_up columns.")
}

clin[, os_time := fifelse(
  !is.na(days_to_death), days_to_death, days_to_last_follow_up
)]

clin[, os_event := fifelse(
  vital_status == "Dead", 1L, 0L
)]

# Keep relevant columns
clin_os <- clin[, .(
  submitter_id,
  os_time,
  os_event,
  vital_status,
  age_at_diagnosis,
  tumor_stage
)]

# Merge OS + ImmuneScore groups
merged <- merge(immune_dt, clin_os, by = "submitter_id", all.x = TRUE)

out_immune_os <- file.path(out_dir, "immune_hot_cold_TNBC_like_TCGA_OS.tsv")
message("Writing immune + OS table to: ", out_immune_os)
fwrite(merged, out_immune_os, sep = "\t")

message("Final rows: ", nrow(merged))
message("Preview:")
print(head(merged))
