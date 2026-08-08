#!/usr/bin/env Rscript

message("=== Phase I – Immune Deconvolution (xCell on TNBC-like TCGA) ===")

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table", repos = "https://cloud.r-project.org")
  }
  library(data.table)
})

# ---- Install Bioc + org.Hs.eg.db if needed ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db", ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ---- Install xCell if needed ----
if (!requireNamespace("xCell", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cloud.r-project.org")
  }
  devtools::install_github("dviraran/xCell")
}
suppressPackageStartupMessages(library(xCell))

# ---------- Paths ----------
expr_path <- file.path("inputs", "expression", "tcga_expr_vst.tsv")
clin_path <- file.path("inputs", "clinical", "tcga_brca_clinical_case12.csv")
out_dir   <- file.path("results", "immune_scores")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---------- 1) Read expression ----------
message("Reading expression matrix from: ", expr_path)
expr_dt <- data.table::fread(expr_path)

# First column = Ensembl IDs with version
gene_ids_ensembl_version <- expr_dt[[1]]
sample_ids_raw <- colnames(expr_dt)[-1]

expr_mat_raw <- as.matrix(expr_dt[, -1])
rownames(expr_mat_raw) <- gene_ids_ensembl_version
colnames(expr_mat_raw) <- sample_ids_raw

message("Expression matrix (raw): ", nrow(expr_mat_raw), " genes x ", ncol(expr_mat_raw), " samples")
message("Example expression sample IDs (raw): ",
        paste(head(sample_ids_raw, 3), collapse = ", "))

# ---------- 2) Convert Ensembl -> SYMBOL ----------
message("Converting Ensembl IDs to HGNC symbols for xCell...")

# strip version: ENSG00000141867.18 -> ENSG00000141867
gene_ids_nover <- sub("\\..*$", "", gene_ids_ensembl_version)

# Build table for mapping
expr_dt_map <- data.table(
  ENSEMBL = gene_ids_nover,
  expr_dt[, -1]
)

# Get mapping from org.Hs.eg.db
map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(expr_dt_map$ENSEMBL),
  keytype = "ENSEMBL",
  columns = "SYMBOL"
)

# Remove rows without SYMBOL
map <- map[!is.na(map$SYMBOL), ]

# Merge mapping into expression
merged_dt <- merge(
  map,
  expr_dt_map,
  by = "ENSEMBL"
)

# Convert to data.table for aggregation
merged_dt <- as.data.table(merged_dt)

# Now aggregate by SYMBOL (some Ensembl IDs share a gene symbol)
expr_cols <- setdiff(colnames(merged_dt), c("ENSEMBL", "SYMBOL"))

agg_dt <- merged_dt[, lapply(.SD, mean), by = SYMBOL, .SDcols = expr_cols]

# Final gene-symbol matrix for xCell
expr_mat <- as.matrix(agg_dt[, ..expr_cols])
rownames(expr_mat) <- agg_dt$SYMBOL
colnames(expr_mat) <- expr_cols  # these are the same original sample IDs

message("Expression matrix (SYMBOL-based): ", nrow(expr_mat), " genes x ", ncol(expr_mat), " samples")

# ---------- 3) Derive case IDs from sample IDs ----------
sample_case_ids <- substr(colnames(expr_mat), 1, 12)
message("Example expression case IDs (12-char): ",
        paste(head(sample_case_ids, 3), collapse = ", "))

# ---------- 4) Read clinical (TCGA case12) ----------
message("Reading TCGA clinical file from: ", clin_path)
clin <- read.csv(clin_path, stringsAsFactors = FALSE, check.names = FALSE)

message("Clinical columns found: ", paste(colnames(clin), collapse = ", "))

possible_cols <- c("submitter_id", "case_id", "patient_id",
                   "barcode", "PATIENT_ID", "Sample_ID", "SAMPLE_ID", "sample_id")

hit <- possible_cols[possible_cols %in% colnames(clin)]

if (length(hit) > 0) {
  sample_col <- hit[1]
  message("Using clinical sample column: ", sample_col)
} else {
  tcga_like_cols <- colnames(clin)[sapply(clin, function(x) {
    any(grepl("^TCGA-", as.character(x)))
  })]
  if (length(tcga_like_cols) > 0) {
    sample_col <- tcga_like_cols[1]
    message("Detected TCGA-like ID column: ", sample_col)
  } else {
    sample_col <- colnames(clin)[1]
    warning("Could not automatically detect a TCGA ID column. ",
            "Falling back to FIRST column: ", sample_col)
  }
}

clin_ids <- clin[[sample_col]]
clin_ids <- gsub("\\.", "-", as.character(clin_ids))

message("Example clinical IDs from column '", sample_col, "': ",
        paste(head(clin_ids, 3), collapse = ", "))

# ---------- 5) Match case IDs ----------
common_cases <- intersect(clin_ids, sample_case_ids)

if (length(common_cases) == 0) {
  stop("No overlap between TCGA clinical IDs (column '", sample_col,
       "') and expression case IDs (first 12 chars of expression column names).\n",
       "Check a few examples manually:\n",
       "  - clinical IDs (head): ", paste(head(clin_ids, 5), collapse = ", "), "\n",
       "  - expression case IDs (head): ", paste(head(sample_case_ids, 5), collapse = ", "))
}

message("Samples with both expr+clin (case-level): ", length(common_cases))

case_to_col <- match(common_cases, sample_case_ids)
expr_sub <- expr_mat[, case_to_col, drop = FALSE]
colnames(expr_sub) <- common_cases

message("Sub-matrix used for xCell: ", nrow(expr_sub), " genes x ",
        ncol(expr_sub), " samples")

# ---------- 6) Run xCell ----------
message("Running xCellAnalysis (this may take a while on first run)...")
xcell_scores <- xCellAnalysis(expr_sub)

out_file <- file.path(out_dir, "xcell_scores_TNBC_like_TCGA.csv")
message("Writing xCell scores to: ", out_file)
write.csv(xcell_scores, out_file, quote = TRUE)

message("Done. Score matrix dimensions: ")
print(dim(xcell_scores))

message("Preview (first few cell types x first few samples):")
print(head(xcell_scores[, 1:min(5, ncol(xcell_scores))]))
