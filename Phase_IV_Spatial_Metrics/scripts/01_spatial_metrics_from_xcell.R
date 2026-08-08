#!/usr/bin/env Rscript

cat("=== Phase IV – Virtual Spatial Metrics from xCell (TNBC-like TCGA) ===\n")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# ---------- Paths ----------
xcell_file <- file.path("..", "Phase_I_Immune_Deconv", "results", "immune_scores",
                        "xcell_scores_TNBC_like_TCGA.csv")

core_file  <- file.path("..", "Phase_II_Immune_Subtypes", "results", "immune_subtypes",
                        "immune_subtypes_TNBC_like_TCGA_core.tsv")

out_dir    <- file.path("results", "spatial_metrics")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Reading xCell scores from:\n  ", xcell_file, "\n")
xcell_dt <- fread(xcell_file)

cat("xCell raw dimensions:", nrow(xcell_dt), "cell types x", ncol(xcell_dt) - 1, "samples\n")

# First column is cell type name
colnames(xcell_dt)[1] <- "cell_type"

# Convert to matrix with cell types as rownames
xcell_mat <- as.matrix(xcell_dt[, -1, with = FALSE])
rownames(xcell_mat) <- xcell_dt$cell_type

# -----------------------------
# Define Immune & Stromal cell sets (xCell names)
# -----------------------------
immune_candidates <- c(
  "CD8+ T-cells", "CD4+ memory T-cells", "CD4+ naive T-cells",
  "Th1 cells", "Th2 cells", "Tregs",
  "NK cells", "NKT",
  "B-cells",
  "Macrophages", "M1 Macrophages", "M2 Macrophages",
  "Monocytes", "Neutrophils",
  "aDC", "cDC", "pDC", "Dendritic cells"
)

stroma_candidates <- c(
  "Fibroblasts",
  "Adipocytes",
  "Endothelial cells",
  "Pericytes",
  "Smooth muscle"
)

immune_present <- intersect(immune_candidates, rownames(xcell_mat))
stroma_present <- intersect(stroma_candidates, rownames(xcell_mat))

cat("Immune cell types present:", paste(immune_present, collapse = ", "), "\n")
cat("Stromal cell types present:", paste(stroma_present, collapse = ", "), "\n")

if (length(immune_present) == 0L) {
  stop("No immune cell types found in xCell matrix with the candidate list.")
}
if (length(stroma_present) == 0L) {
  stop("No stromal cell types found in xCell matrix with the candidate list.")
}

# Compute composite scores per sample
immune_scores <- colMeans(xcell_mat[immune_present, , drop = FALSE])
stroma_scores <- colMeans(xcell_mat[stroma_present, , drop = FALSE])

spatial_dt <- data.table(
  submitter_id   = colnames(xcell_mat),
  Immune_composite = as.numeric(immune_scores[colnames(xcell_mat)]),
  Stroma_composite = as.numeric(stroma_scores[colnames(xcell_mat)])
)

# Avoid division by zero
eps <- 1e-6

spatial_dt[, Immune_Stroma_ratio :=
             Immune_composite / (Immune_composite + Stroma_composite + eps)]

spatial_dt[, Immune_Exclusion_index :=
             Stroma_composite / (Immune_composite + Stroma_composite + eps)]

cat("Spatial metrics summary (head):\n")
print(head(spatial_dt))

# -----------------------------
# Merge with immune subtype + OS table
# -----------------------------
cat("Reading core immune subtype table from:\n  ", core_file, "\n")
core_dt <- fread(core_file)

cat("Core table dimensions:", nrow(core_dt), "x", ncol(core_dt), "\n")
cat("Core columns:", paste(colnames(core_dt), collapse = ", "), "\n")

merged_dt <- merge(core_dt, spatial_dt, by = "submitter_id", all.x = TRUE)

cat("Merged spatial + clinical dimensions:", nrow(merged_dt), "x", ncol(merged_dt), "\n")

# Create tertiles for Immune_Stroma_ratio
q <- quantile(merged_dt$Immune_Stroma_ratio, probs = c(1/3, 2/3), na.rm = TRUE)
cat("Immune_Stroma_ratio tertiles:", paste(round(q, 4), collapse = " / "), "\n")

merged_dt[, Immune_Stroma_tertile :=
            fifelse(Immune_Stroma_ratio <= q[1], "Low",
              fifelse(Immune_Stroma_ratio <= q[2], "Mid", "High"))]

merged_dt$Immune_Stroma_tertile <- factor(
  merged_dt$Immune_Stroma_tertile,
  levels = c("Low", "Mid", "High")
)

# Overall summary
cat("Immune_Stroma_tertile counts:\n")
print(table(merged_dt$Immune_Stroma_tertile, useNA = "ifany"))

# -----------------------------
# Save tables
# -----------------------------
raw_out   <- file.path(out_dir, "spatial_metrics_raw_TNBC_like_TCGA.tsv")
merged_out <- file.path(out_dir, "spatial_metrics_with_clin_TNBC_like_TCGA.tsv")

fwrite(spatial_dt, raw_out, sep = "\t")
fwrite(merged_dt, merged_out, sep = "\t")

cat("Saved raw spatial metrics to:\n  ", raw_out, "\n")
cat("Saved spatial + clinical merged table to:\n  ", merged_out, "\n")

# -----------------------------
# Simple QC plots
# -----------------------------
p1 <- ggplot(merged_dt, aes(x = Immune_Stroma_ratio)) +
  geom_histogram(bins = 30) +
  theme_bw() +
  ggtitle("Distribution of Immune:Stroma ratio (TNBC-like TCGA)")

ggsave(file.path(out_dir, "hist_Immune_Stroma_ratio_TNBC_like_TCGA.png"),
       p1, width = 6, height = 4, dpi = 300)

p2 <- ggplot(merged_dt, aes(x = immune_subtype, y = Immune_Stroma_ratio,
                            fill = immune_subtype)) +
  geom_boxplot(outlier.size = 0.8) +
  theme_bw() +
  xlab("Immune subtype") +
  ylab("Immune:Stroma ratio") +
  ggtitle("Immune:Stroma ratio by immune subtype (TNBC-like TCGA)")

ggsave(file.path(out_dir, "boxplot_Immune_Stroma_ratio_by_subtype_TNBC_like_TCGA.png"),
       p2, width = 6, height = 4, dpi = 300)

cat("Saved QC plots to:", out_dir, "\n")
cat("=== DONE Phase IV – Spatial metrics build ===\n")