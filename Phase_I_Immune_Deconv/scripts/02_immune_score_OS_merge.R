#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

message("=== Phase I – ImmuneScore + OS merge (TNBC-like TCGA) ===")

xcell_path <- "results/immune_scores/xcell_scores_TNBC_like_TCGA.csv"
clin_path  <- "inputs/clinical/tcga_brca_clinical_case12.csv"
out_sample <- "results/immune_scores/xcell_scores_TNBC_like_TCGA_samples.tsv"
out_immune <- "results/immune_scores/immune_hot_cold_TNBC_like_TCGA_OS.tsv"

# ---------------------------
# 1) Read xCell scores (cell_types x samples)
# ---------------------------
message("Reading xCell scores from: ", xcell_path)
xcell <- as.matrix(read.csv(xcell_path, row.names = 1, check.names = FALSE))
message("xCell matrix: ", paste(dim(xcell), collapse = " x "))

# transpose: samples x cell_types
xcell_t <- as.data.table(t(xcell), keep.rownames = "submitter_id")

message("Writing sample-wise xCell table to: ", out_sample)
dir.create(dirname(out_sample), recursive = TRUE, showWarnings = FALSE)
fwrite(xcell_t, out_sample, sep = "\t")

# ---------------------------
# 2) Use xCell ImmuneScore row
# ---------------------------
if (!"ImmuneScore" %in% rownames(xcell)) {
  stop("ImmuneScore row not found in xCell output.")
}

immune_vec <- xcell["ImmuneScore", ]
immune_dt <- data.table(
  submitter_id = names(immune_vec),
  ImmuneScore  = as.numeric(immune_vec)
)

message("ImmuneScore summary:")
print(summary(immune_dt$ImmuneScore))

# Tertiles: Cold / Intermediate / Hot
q1 <- quantile(immune_dt$ImmuneScore, 1/3, na.rm = TRUE)
q2 <- quantile(immune_dt$ImmuneScore, 2/3, na.rm = TRUE)

immune_dt[, immune_group := fifelse(
  ImmuneScore <= q1, "Cold",
  fifelse(ImmuneScore >= q2, "Hot", "Intermediate")
)]

message("Group counts:\n")
print(table(immune_dt$immune_group))

# ---------------------------
# 3) Merge with OS clinical
# ---------------------------
message("Reading TCGA clinical file from: ", clin_path)
clin <- fread(clin_path)

needed_cols <- c("submitter_id", "os_time", "os_event")
missing <- setdiff(needed_cols, names(clin))
if (length(missing) > 0) {
  stop("Clinical file missing columns: ", paste(missing, collapse = ", "))
}

merged <- merge(immune_dt, clin[, ..needed_cols], by = "submitter_id", all.x = TRUE)

# Optional: basic type cleaning
merged[, os_time  := as.numeric(os_time)]
merged[, os_event := as.integer(os_event)]

message("Writing immune + OS table to: ", out_immune)
dir.create(dirname(out_immune), recursive = TRUE, showWarnings = FALSE)
fwrite(merged, out_immune, sep = "\t")

message("Final rows: ", nrow(merged))
message("Preview:")
print(head(merged, 6))

message("=== DONE ImmuneScore OS merge ===")
