library(data.table)

input_file <- "../Phase_IV_Spatial_Metrics/results/spatial_metrics/spatial_metrics_with_clin_TNBC_like_TCGA.tsv"
out_file   <- "results/drug_penetration/DPI_TNBC_like_TCGA.tsv"

cat("=== Phase V – Compute Drug Penetration Index (DPI) ===\n")

dt <- fread(input_file)
cat("Loaded rows:", nrow(dt), "\n")

# ---- DPI Formula ----
dt[, DPI := 0.6 * (1 - Immune_Exclusion_index) +
           0.4 * (1 - Stroma_composite)]

# ---- Tertiles ----
dt[, DPI_tertile := cut(DPI,
                        breaks = quantile(DPI, probs = c(0,0.33,0.66,1)),
                        include.lowest = TRUE,
                        labels = c("Low","Mid","High"))]

fwrite(dt, out_file, sep="\t")
cat("Saved:", out_file, "\n")

summary(dt$DPI)
table(dt$DPI_tertile)