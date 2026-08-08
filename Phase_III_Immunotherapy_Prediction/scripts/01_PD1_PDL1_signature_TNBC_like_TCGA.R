#!/usr/bin/env Rscript

cat("=== Phase III – PD1/PDL1 Signature (TNBC-like TCGA) ===\n")

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(dplyr)
})

## --- Paths (relative to Phase_III root) ---
expr_file      <- "inputs/expression/tcga_expr_vst.tsv"
immune_os_file <- "../Phase_I_Immune_Deconv/results/immune_scores/immune_hot_cold_TNBC_like_TCGA_OS.tsv"
out_dir        <- "results/immunotherapy"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Reading expression matrix from:", expr_file, "\n")
expr_dt <- fread(expr_file)
cat("Expression matrix (raw):", nrow(expr_dt), "genes x", ncol(expr_dt) - 1, "samples\n")

gene_ids <- expr_dt[[1]]
expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
rownames(expr_mat) <- gene_ids
sample_ids_raw <- colnames(expr_mat)
case_ids       <- substr(sample_ids_raw, 1, 12)
cat("Example expression sample IDs (raw):", paste(head(sample_ids_raw, 3), collapse = ", "), "\n")

## --- Map Ensembl → SYMBOL ---
cat("Converting Ensembl IDs to HGNC symbols...\n")
ens_base <- sub("\\..*", "", rownames(expr_mat))
map_df <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(ens_base),
  keytype = "ENSEMBL",
  columns = "SYMBOL"
)
map_dt <- as.data.table(map_df)[!is.na(SYMBOL)]

merged_dt <- data.table(ENSEMBL = ens_base, expr_mat)
merged_dt <- merge(map_dt, merged_dt, by = "ENSEMBL")

expr_cols <- setdiff(names(merged_dt), c("ENSEMBL", "SYMBOL"))
symbol_dt <- merged_dt[, lapply(.SD, mean), by = SYMBOL, .SDcols = expr_cols]

symbol_mat <- as.matrix(symbol_dt[, -1, with = FALSE])
rownames(symbol_mat) <- symbol_dt$SYMBOL
colnames(symbol_mat) <- expr_cols

cat("Expression matrix (SYMBOL-based):", nrow(symbol_mat), "genes x", ncol(symbol_mat), "samples\n")

## --- Define PD1/PDL1 / ICB signature genes ---
sig_genes <- c(
  "CD274",   # PD-L1
  "PDCD1",   # PD-1
  "CTLA4",
  "LAG3",
  "TIGIT",
  "CXCL9",
  "CXCL10",
  "CXCL13",
  "IFNG",
  "GZMB"
)

present_genes <- intersect(sig_genes, rownames(symbol_mat))
cat("Signature genes present in matrix:", paste(present_genes, collapse = ", "), "\n")

if (length(present_genes) < 3) {
  stop("Too few signature genes found in expression matrix. Found only: ",
       paste(present_genes, collapse = ", "))
}

sig_mat <- symbol_mat[present_genes, , drop = FALSE]

## --- Row-wise z-score then average = signature ---
cat("Computing per-gene z-scores and overall PD1/PDL1 signature...\n")
gene_z <- t(scale(t(sig_mat)))
signature_scores <- colMeans(gene_z, na.rm = TRUE)

sig_dt <- data.table(
  sample_id = colnames(symbol_mat),
  case_id   = substr(colnames(symbol_mat), 1, 12),
  PD1_PDL1_signature = as.numeric(signature_scores)
)

## --- Merge with TNBC-like OS table (Phase I) ---
cat("Reading TNBC-like immune + OS file from:", immune_os_file, "\n")
immune_os <- fread(immune_os_file)

if (!"submitter_id" %in% names(immune_os)) {
  stop("submitter_id column not found in ", immune_os_file)
}

cat("TNBC-like immune+OS rows:", nrow(immune_os), "\n")

merged <- merge(
  immune_os,
  unique(sig_dt[, .(case_id, PD1_PDL1_signature)]),
  by.x = "submitter_id",
  by.y = "case_id",
  all.x = TRUE
)

cat("Merged TNBC-like rows with PD1/PDL1 signature:", nrow(merged), "\n")
cat("Non-NA signature rows:", sum(!is.na(merged$PD1_PDL1_signature)), "\n")

## --- Save table ---
out_table <- file.path(out_dir, "PD1_PDL1_signature_TNBC_like_TCGA.tsv")
fwrite(merged, out_table, sep = "\t")
cat("Wrote PD1/PDL1 signature table to:", out_table, "\n")

## --- Survival analysis: median split high vs low ---
surv_dt <- merged %>%
  filter(!is.na(os_time), !is.na(os_event), !is.na(PD1_PDL1_signature))

cat("Survival subset rows:", nrow(surv_dt), "\n")

if (nrow(surv_dt) > 0) {
  med_sig <- median(surv_dt$PD1_PDL1_signature, na.rm = TRUE)
  surv_dt <- surv_dt %>%
    mutate(
      PD1_PDL1_group = ifelse(PD1_PDL1_signature >= med_sig, "High", "Low")
    )

  fit <- survfit(Surv(os_time, os_event) ~ PD1_PDL1_group, data = surv_dt)
  p_kw <- surv_pvalue(fit, data = surv_dt)$pval

  cat("KM p-value (High vs Low PD1/PDL1 signature):", p_kw, "\n")

  p_km <- ggsurvplot(
    fit,
    data           = surv_dt,
    pval           = TRUE,
    risk.table     = TRUE,
    ggtheme        = theme_bw(),
    legend.title   = "PD1/PDL1 signature",
    legend.labs    = c("High", "Low")
  )

  km_file <- file.path(out_dir, "KM_PD1_PDL1_signature_high_vs_low_TNBC_like_TCGA.png")
  ggsave(km_file, p_km$plot, width = 6, height = 5, dpi = 300)
  cat("Saved KM plot to:", km_file, "\n")

  ## also save risk table as separate PNG to avoid ggplot2+table clash
  rt_file <- file.path(out_dir, "KM_PD1_PDL1_signature_high_vs_low_TNBC_like_TCGA_risktable.png")
  ggsave(rt_file, p_km$table, width = 6, height = 3, dpi = 300)
  cat("Saved KM risk table plot to:", rt_file, "\n")
} else {
  cat("Not enough rows for survival analysis; skipping KM/Cox.\n")
}

cat("=== DONE PD1/PDL1 Signature Phase ===\n")
