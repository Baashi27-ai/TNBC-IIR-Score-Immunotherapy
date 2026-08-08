#!/usr/bin/env Rscript

cat("=== Phase VI – Download TCGA BRCA MAF (Mutect2) ===\n")

suppressPackageStartupMessages({
  library(data.table)
  library(TCGAbiolinks)
})

# Check function exists
if (!"GDCquery_Maf" %in% ls("package:TCGAbiolinks")) {
  stop("ERROR: TCGAbiolinks is loaded but function GDCquery_Maf() is not found.\n",
       "Your version may be older. Run this in R:\n",
       "BiocManager::install('TCGAbiolinks', force=TRUE)\n")
}

out_dir  <- "../inputs"
maf_file <- file.path(out_dir, "TCGA_BRCA_mutect2_MAF.tsv.gz")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

if (file.exists(maf_file)) {
  cat("MAF already exists at:\n  ", maf_file, "\nSkipping download.\n")
  quit(save="no", status = 0)
}

cat("Downloading Mutect2 MAF from GDC… this may take a few minutes…\n")

# --- Actual MAF download using Mutect2 ---
maf <- TCGAbiolinks::GDCquery_Maf(
  tumor = "BRCA",
  pipelines = "mutect2"
)

cat("Downloaded rows:", nrow(maf), "\n")

cat("Saving to:", maf_file, "\n")
fwrite(maf, maf_file, sep="\t")

cat("=== DONE: BRCA Mutect2 MAF downloaded ===\n")
