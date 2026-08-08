#!/usr/bin/env Rscript

cat("=== Phase V – Therapy Accessibility Mapping (TNBC-like TCGA) ===\n")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Inputs
spatial_file <- "../Phase_IV_Spatial_Metrics/results/spatial_metrics/spatial_metrics_with_clin_TNBC_like_TCGA.tsv"
dpi_file     <- "results/drug_penetration/DPI_TNBC_like_TCGA.tsv"

# Outputs
out_table <- "results/drug_penetration/therapy_accessibility_TNBC_like_TCGA.tsv"
out_bar   <- "results/drug_penetration/therapy_accessibility_counts_TNBC_like_TCGA.png"

## --- Load inputs ------------------------------------------------------------
if (!file.exists(spatial_file)) {
  stop("Spatial metrics file not found: ", spatial_file)
}
if (!file.exists(dpi_file)) {
  stop("DPI file not found: ", dpi_file)
}

spat <- fread(spatial_file)
dpi  <- fread(dpi_file)

cat("Spatial table:", nrow(spat), "rows x", ncol(spat), "cols\n")
cat("DPI table    :", nrow(dpi), "rows x", ncol(dpi), "cols\n")

## --- Harmonise column names -------------------------------------------------
# Fix capitalization if needed
if ("Immune_exclusion_index" %in% names(spat)) {
  setnames(spat, "Immune_exclusion_index", "Immune_Exclusion_index")
}
if (!"Immune_Exclusion_index" %in% names(spat)) {
  stop("Column 'Immune_Exclusion_index' not found in spatial metrics table.")
}

# Which DPI continuous column?
dpi_cont_candidates <- intersect(c("DPI_index", "DPI_score", "DPI", "DPI_composite"), names(dpi))
if (length(dpi_cont_candidates) == 0) {
  stop("Could not find a continuous DPI column in DPI file.")
}
dpi_cont_col <- dpi_cont_candidates[1]
cat("Using DPI continuous column:", dpi_cont_col, "\n")

# Ensure DPI_tertile exists
if (!"DPI_tertile" %in% names(dpi)) {
  cat("DPI_tertile not found; creating tertiles from", dpi_cont_col, "\n")
  q <- quantile(dpi[[dpi_cont_col]], probs = c(1/3, 2/3), na.rm = TRUE)
  dpi[, DPI_tertile := cut(
    get(dpi_cont_col),
    breaks = c(-Inf, q[1], q[2], Inf),
    labels = c("Low", "Mid", "High"),
    include.lowest = TRUE
  )]
}

## --- Merge spatial + DPI on submitter_id -----------------------------------
common_id <- intersect(spat$submitter_id, dpi$submitter_id)
cat("Common submitter_id rows:", length(common_id), "\n")

dpi_sub <- dpi[, .(submitter_id,
                   DPI_index   = get(dpi_cont_col),
                   DPI_tertile = DPI_tertile)]

merged <- merge(spat, dpi_sub, by = "submitter_id", all.x = TRUE)

# Ensure numeric
merged[, Immune_Exclusion_index := as.numeric(Immune_Exclusion_index)]
merged[, DPI_index              := as.numeric(DPI_index)]

## --- Exclusion tertiles -----------------------------------------------------
if (!"Exclusion_tertile" %in% names(merged)) {
  cat("Creating Exclusion_tertile from Immune_Exclusion_index\n")
  q_ex <- quantile(merged$Immune_Exclusion_index, probs = c(1/3, 2/3), na.rm = TRUE)
  merged[, Exclusion_tertile := cut(
    Immune_Exclusion_index,
    breaks = c(-Inf, q_ex[1], q_ex[2], Inf),
    labels = c("Low", "Mid", "High"),
    include.lowest = TRUE
  )]
}

## --- Therapy accessibility groups ------------------------------------------
# Heuristic:
#   High_access:  DPI High & Exclusion Low
#   Poor_access:  DPI Low  & Exclusion High
#   Medium:       everything else

merged[, therapy_access_group := "Medium"]
merged[DPI_tertile == "High" & Exclusion_tertile == "Low",  therapy_access_group := "High_access"]
merged[DPI_tertile == "Low"  & Exclusion_tertile == "High", therapy_access_group := "Poor_access"]

merged[, therapy_access_group := factor(
  therapy_access_group,
  levels = c("High_access", "Medium", "Poor_access")
)]

## --- Save table -------------------------------------------------------------
dir.create(dirname(out_table), showWarnings = FALSE, recursive = TRUE)
fwrite(merged, out_table, sep = "\t")
cat("Saved therapy accessibility table to:", out_table, "\n")

## --- Simple barplot of group counts ----------------------------------------
tab <- merged[!is.na(therapy_access_group), .N, by = therapy_access_group]
cat("Therapy access group counts:\n")
print(tab)

p_bar <- ggplot(tab, aes(x = therapy_access_group, y = N, fill = therapy_access_group)) +
  geom_col() +
  theme_minimal() +
  labs(
    x = "Therapy accessibility group",
    y = "Number of patients",
    title = "Estimated therapy accessibility (DPI + Immune exclusion)"
  ) +
  theme(legend.position = "none")

ggsave(out_bar, p_bar, width = 5, height = 4, dpi = 300)
cat("Saved barplot to:", out_bar, "\n")

cat("=== DONE therapy accessibility mapping ===\n")