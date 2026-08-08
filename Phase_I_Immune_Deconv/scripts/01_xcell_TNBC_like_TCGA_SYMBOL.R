#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(xCell)
})

message("=== Phase I – Immune Deconvolution (xCell on TNBC-like TCGA) ===")

expr_path  <- "inputs/expression/tcga_expr_vst.tsv"
clin_path  <- "inputs/clinical/tcga_brca_clinical_case12.csv"
out_xcell  <- "results/immune_scores/xcell_scores_TNBC_like_TCGA.csv"

# ---------------------------
# 1) Load expression (ENSEMBL)
# ---------------------------
message("Reading expression matrix from: ", expr_path)
expr_dt <- fread(expr_path)
message("Expression matrix (raw): ", paste(dim(expr_dt), collapse = " x "))

gene_col <- names(expr_dt)[1]
expr_mat_raw <- as.matrix(expr_dt[, -1])
rownames(expr_mat_raw) <- expr_dt[[gene_col]]

message("Example expression sample IDs (raw): ",
        paste(colnames(expr_mat_raw)[1:3], collapse = ", "))

# ---------------------------
# 2) Map ENSEMBL → SYMBOL
# ---------------------------
message("Converting Ensembl IDs to HGNC symbols for xCell...")

ens_ids <- sub("\\.\\d+$", "", rownames(expr_mat_raw))  # strip .version
map_df <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = unique(ens_ids),
  keytype  = "ENSEMBL",
  columns  = c("SYMBOL")
)

map_dt <- as.data.table(map_df)
setnames(map_dt, c("ENSEMBL", "SYMBOL"))

# Merge mapping with expression
expr_dt_long <- data.table(ENSEMBL = ens_ids)
expr_dt_long[, row_id := .I]
expr_dt_long <- cbind(expr_dt_long, as.data.table(expr_mat_raw))

merged_dt <- merge(map_dt[!is.na(SYMBOL)], expr_dt_long, by = "ENSEMBL")

expr_cols <- setdiff(names(merged_dt), c("ENSEMBL", "SYMBOL", "row_id"))

# Average rows with same SYMBOL
expr_symbol_dt <- merged_dt[, lapply(.SD, mean), by = SYMBOL, .SDcols = expr_cols]

expr_symbol_mat <- as.matrix(expr_symbol_dt[, ..expr_cols])
rownames(expr_symbol_mat) <- expr_symbol_dt$SYMBOL

message("Expression matrix (SYMBOL-based): ",
        paste(dim(expr_symbol_mat), collapse = " x "))

# ---------------------------
# 3) Align with clinical IDs
# ---------------------------
message("Example expression case IDs (12-char): ",
        paste(substr(colnames(expr_symbol_mat)[1:3], 1, 12), collapse = ", "))

message("Reading TCGA clinical file from: ", clin_path)
clin <- fread(clin_path)
message("Clinical columns found: ", paste(names(clin), collapse = ", "))

# We expect submitter_id + os_time + os_event from OSbuild
if (!"submitter_id" %in% names(clin)) {
  stop("Clinical file must contain 'submitter_id' column after OSbuild.")
}

clin[, case_id := substr(submitter_id, 1, 12)]

expr_case_ids <- substr(colnames(expr_symbol_mat), 1, 12)
clin_case_ids <- clin$case_id

common_ids <- intersect(expr_case_ids, clin_case_ids)
message("Samples with both expr+clin (case-level): ", length(common_ids))

if (length(common_ids) < 10) {
  stop("Too few overlapping samples between expression and clinical; check IDs.")
}

keep_cols <- which(expr_case_ids %in% common_ids)
expr_symbol_sub <- expr_symbol_mat[, keep_cols, drop = FALSE]

message("Sub-matrix used for xCell: ",
        paste(dim(expr_symbol_sub), collapse = " x "))

# Rename cols to 12-char case IDs (xCell doesn’t care)
colnames(expr_symbol_sub) <- substr(colnames(expr_symbol_sub), 1, 12)

# ---------------------------
# 4) Run xCell
# ---------------------------
message("Running xCellAnalysis (this may take a while on first run)...")
xcell_scores <- xCell::xCellAnalysis(expr_symbol_sub)

# xCell returns: cell_types x samples
message("Writing xCell scores to: ", out_xcell)
dir.create(dirname(out_xcell), recursive = TRUE, showWarnings = FALSE)
write.csv(xcell_scores, out_xcell, quote = FALSE)

message("Done. Score matrix dimensions:")
print(dim(xcell_scores))

message("Preview (first few cell types x first few samples):")
print(xcell_scores[1:6, 1:4])

message("=== DONE xCell TNBC-like TCGA ===")
