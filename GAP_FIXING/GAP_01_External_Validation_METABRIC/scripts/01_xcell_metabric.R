#!/usr/bin/env Rscript
################################################################################
# Script: 01_xcell_metabric.R
# Purpose: Run xCell deconvolution on METABRIC TNBC expression data
# Input:  data/METABRIC_TNBC_expression.tsv
# Output: results/xcell/xcell_scores_metabric.tsv
################################################################################

cat("\n")
cat("============================================================\n")
cat(" xCELL DECONVOLUTION - METABRIC TNBC (R version)\n")
cat("============================================================\n")
cat("Started at:", Sys.time(), "\n\n")

# =============================================================================
# 1. SETUP
# =============================================================================

# Set working directory
setwd("D:/Baashi/TNBC_project/Immune_Spatial_Immunotherapy/GAP_FIXING/GAP_01_External_Validation_METABRIC")

# Create results directories
dir.create("results/xcell", recursive = TRUE, showWarnings = FALSE)
dir.create("data", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 2. INSTALL/LOAD PACKAGES
# =============================================================================

cat("Checking required packages...\n")

# Check if xCell is installed
if (!requireNamespace("xCell", quietly = TRUE)) {
    cat("  xCell package not found. Installing from GitHub...\n")
    
    # Install devtools if needed
    if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("devtools", repos = "https://cran.rstudio.com")
    }
    
    # Install xCell from GitHub
    devtools::install_github("dviraran/xCell")
}

# Load packages
library(xCell)
library(data.table)

cat("  ✓ Packages loaded\n\n")

# =============================================================================
# 3. LOAD DATA
# =============================================================================

cat("Loading METABRIC expression data...\n")

expr_file <- "../../../M9_external_validation/inputs/METABRIC_TNBC_expression.tsv"

if (!file.exists(expr_file)) {
    stop("Expression file not found: ", expr_file)
}

# Load expression
expr_df <- fread(expr_file, header = TRUE, data.table = FALSE)

cat("  Expression matrix:", nrow(expr_df), "genes ×", ncol(expr_df) - 1, "samples\n")

# Convert to matrix (genes as rownames)
rownames(expr_df) <- expr_df$Hugo_Symbol
expr_matrix <- as.matrix(expr_df[, -1])

# Check for duplicate genes
if (any(duplicated(rownames(expr_matrix)))) {
    cat("  Warning: Duplicate gene symbols found. Keeping first occurrence.\n")
    expr_matrix <- expr_matrix[!duplicated(rownames(expr_matrix)), ]
}

cat("  Expression matrix ready:", nrow(expr_matrix), "genes ×", ncol(expr_matrix), "samples\n\n")

# =============================================================================
# 4. RUN xCELL
# =============================================================================

cat("Running xCell deconvolution...\n")
cat("  This will take 5-10 minutes...\n")

# Run xCell with default settings
xcell_result <- xCellAnalysis(
    expr_matrix,
    genes = NULL,
    cell.types = NULL,
    alpha = 0.5,
    spill = FALSE,        # Disable spillover for speed
    parallel = FALSE,     # Disable parallel for stability
    save.raw = FALSE,
    signatures = NULL     # Use default xCell signatures
)

cat("  ✓ xCell deconvolution complete\n\n")

# =============================================================================
# 5. FORMAT RESULTS
# =============================================================================

cat("Formatting results...\n")

# Transpose to get samples as rows, cell types as columns
xcell_df <- as.data.frame(t(xcell_result))

# Add sample IDs
xcell_df$sample_id <- rownames(xcell_df)

# Reorder columns to have sample_id first
xcell_df <- xcell_df[, c(ncol(xcell_df), 1:(ncol(xcell_df)-1))]

# =============================================================================
# 6. SAVE RESULTS
# =============================================================================

output_file <- "results/xcell/xcell_scores_metabric.tsv"

fwrite(xcell_df, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("  ✓ xCell scores saved to:", output_file, "\n")
cat("  Dimensions:", nrow(xcell_df), "samples ×", ncol(xcell_df) - 1, "cell types\n\n")

# =============================================================================
# 7. SUMMARY
# =============================================================================

cat("============================================================\n")
cat(" SUMMARY\n")
cat("============================================================\n\n")

# Get immune cell types (based on xCell documentation)
immune_types <- c(
    "B cells", "CD4+ T-cells", "CD8+ T-cells", "T cells", "Tregs",
    "NK cells", "Macrophages M1", "Macrophages M2", "Monocytes",
    "Myeloid dendritic cells", "Plasmacytoid dendritic cells"
)

available_immune <- intersect(immune_types, colnames(xcell_df))

if (length(available_immune) > 0) {
    cat("Immune cell type scores (mean ± sd):\n")
    for (ct in available_immune) {
        cat(sprintf("  %s: %.4f ± %.4f\n", 
            ct, 
            mean(xcell_df[[ct]], na.rm = TRUE),
            sd(xcell_df[[ct]], na.rm = TRUE)
        ))
    }
}

cat("\nCompleted at:", Sys.time(), "\n")
cat("\n✅ xCELL DECONVOLUTION COMPLETE\n")
cat("   Next: 02_compute_immune_scores.R\n\n")