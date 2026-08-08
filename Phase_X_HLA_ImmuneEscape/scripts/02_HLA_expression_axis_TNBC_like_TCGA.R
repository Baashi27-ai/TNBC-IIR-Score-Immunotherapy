#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

cat("=== Phase X – MHC-I Expression Axis (TNBC-like TCGA) ===\n")

## Always work from project root
setwd("/mnt/c/TNBC_project")

## 1) Files
expr_file     <- "Immune_Spatial_Immunotherapy/Phase_I_Immune_Deconv/inputs/expression/tcga_expr_vst.tsv"
immune_os_file <- "Immune_Spatial_Immunotherapy/Phase_I_Immune_Deconv/results/immune_scores/immune_hot_cold_TNBC_like_TCGA_OS.tsv"

out_dir  <- "Immune_Spatial_Immunotherapy/Phase_X_HLA_ImmuneEscape/results/hla"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "HLA_expression_axis_TNBC_like_TCGA.tsv")

cat("Reading expression matrix from:\n  ", expr_file, "\n")
expr <- fread(expr_file)

cat("Expression dimensions (raw):", nrow(expr), "genes x", ncol(expr) - 1, "samples\n")

gene_col <- names(expr)[1]

## 2) Clean Ensembl IDs
expr[, ENSEMBL := sub("\\..*$", "", get(gene_col))]

## 3) Map ENSEMBL → SYMBOL
ens_ids <- unique(expr$ENSEMBL)

cat("Mapping", length(ens_ids), "Ensembl IDs to HGNC symbols...\n")
map_dt <- as.data.table(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys     = ens_ids,
    keytype  = "ENSEMBL",
    columns  = "SYMBOL"
  )
)

map_dt <- map_dt[!is.na(SYMBOL)]
cat("Mapped to", length(unique(map_dt$SYMBOL)), "unique symbols.\n")

## Merge expression with SYMBOL
expr_long <- merge(
  expr[, c("ENSEMBL", names(expr)[names(expr) != gene_col]), with = FALSE],
  map_dt,
  by = "ENSEMBL",
  all.x = TRUE
)

expr_long <- expr_long[!is.na(SYMBOL)]

## 4) Aggregate by SYMBOL (mean across Ensembl duplicates)
sample_cols <- setdiff(names(expr_long), c("ENSEMBL", "SYMBOL"))
cat("Aggregating expression by SYMBOL...\n")

expr_sym <- expr_long[
  ,
  lapply(.SD, mean, na.rm = TRUE),
  by = SYMBOL,
  .SDcols = sample_cols
]

cat("SYMBOL-based expression:", nrow(expr_sym), "genes x", ncol(expr_sym) - 1, "samples\n")

## 5) Select MHC-I pathway genes
mhc_genes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "B2M", "TAP1", "TAP2")

found_genes <- intersect(mhc_genes, expr_sym$SYMBOL)
cat("MHC-I genes requested:", paste(mhc_genes, collapse = ", "), "\n")
cat("MHC-I genes found in matrix:", paste(found_genes, collapse = ", "), "\n")

if (length(found_genes) == 0L) {
  stop("No requested MHC-I genes found in expression matrix.")
}

mhc_expr <- expr_sym[SYMBOL %in% found_genes]

## 6) Build matrix (genes × samples) and z-score by gene
mat <- as.matrix(mhc_expr[, ..sample_cols])
rownames(mat) <- mhc_expr$SYMBOL

cat("MHC-I expression matrix:", nrow(mat), "genes x", ncol(mat), "samples\n")

## Z-score each gene (row-wise)
mat_z <- t(scale(t(mat)))
colnames(mat_z) <- sample_cols

## 7) Per-sample MHC-I axis = mean z-score across genes
mhc_axis <- colMeans(mat_z, na.rm = TRUE)

axis_dt <- data.table(
  sample_raw   = names(mhc_axis),
  MHC_I_axis   = as.numeric(mhc_axis)
)

## Derive 12-char submitter_id from TCGA barcodes
axis_dt[, submitter_id := substr(sample_raw, 1L, 12L)]

## Aggregate if multiple samples per patient
axis_dt <- axis_dt[
  ,
  .(MHC_I_axis = mean(MHC_I_axis, na.rm = TRUE)),
  by = submitter_id
]

cat("MHC-I axis rows after patient aggregation:", nrow(axis_dt), "\n")

## 8) Restrict to TNBC-like immune OS cohort
cat("Reading immune OS cohort from:\n  ", immune_os_file, "\n")
imm_os <- fread(immune_os_file)

if (!"submitter_id" %in% names(imm_os)) {
  stop("immune OS file does not contain 'submitter_id' column.")
}

tnbc_ids <- unique(imm_os$submitter_id)
axis_dt <- axis_dt[submitter_id %in% tnbc_ids]

cat("MHC-I axis rows after TNBC-like filter:", nrow(axis_dt), "\n")

## 9) Add tertiles
valid <- axis_dt[!is.na(MHC_I_axis)]
if (nrow(valid) > 0L) {
  qs <- quantile(valid$MHC_I_axis, probs = c(1/3, 2/3), na.rm = TRUE)
  axis_dt[
    ,
    MHC_I_axis_tertile := fifelse(
      is.na(MHC_I_axis), NA_character_,
      fifelse(
        MHC_I_axis <= qs[1], "Low",
        fifelse(MHC_I_axis <= qs[2], "Mid", "High")
      )
    )
  ]
} else {
  axis_dt[, MHC_I_axis_tertile := NA_character_]
}

## 10) Save
fwrite(axis_dt, out_file, sep = "\t")

cat("Saved MHC-I axis table to:\n  ", out_file, "\n")
cat("Preview:\n")
print(head(axis_dt))

cat("=== DONE MHC-I Expression Axis ===\n")
