#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(survminer)
  library(ggplot2)
})

cat("=== Phase II – Immune Subtype Summary & Survival (TNBC-like TCGA) ===\n")

## ------------------------------------------------------------------
## 1. Paths
## ------------------------------------------------------------------
core_file <- file.path(
  "Immune_Spatial_Immunotherapy",
  "Phase_II_Immune_Subtypes",
  "results",
  "immune_subtypes",
  "immune_subtypes_TNBC_like_TCGA_core.tsv"
)

out_dir <- file.path(
  "Immune_Spatial_Immunotherapy",
  "Phase_II_Immune_Subtypes",
  "results",
  "immune_subtypes"
)

if (!file.exists(core_file)) {
  stop("Core immune subtype file not found: ", core_file)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------------
## 2. Load core table
## ------------------------------------------------------------------
cat("Reading core immune subtype table from:", core_file, "\n")
dt <- fread(core_file)
cat("Core table dimensions:", paste(dim(dt), collapse = " x "), "\n")
cat("Columns:", paste(names(dt), collapse = ", "), "\n")

## Ensure immune_subtype factor with canonical order
if (!"immune_subtype" %in% names(dt)) {
  stop("Column 'immune_subtype' not found in core table.")
}

dt[, immune_subtype := factor(
  immune_subtype,
  levels = c("BLIS", "IM", "BLIA", "LAR")
)]

## ------------------------------------------------------------------
## 3. Boxplot: ImmuneScore by immune_subtype
## ------------------------------------------------------------------
boxplot_out <- file.path(
  out_dir,
  "boxplot_ImmuneScore_by_immune_subtype_TNBC_like_TCGA.png"
)

p_box <- ggplot(
  dt[!is.na(immune_subtype)],
  aes(x = immune_subtype, y = ImmuneScore, fill = immune_subtype)
) +
  geom_boxplot(outlier.alpha = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
  labs(
    x = "Immune subtype",
    y = "ImmuneScore (xCell)",
    title = "ImmuneScore distribution by immune subtype (TNBC-like TCGA)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

ggsave(boxplot_out, p_box, width = 7, height = 5, dpi = 150)
cat("Saved boxplot to:", boxplot_out, "\n")

## ------------------------------------------------------------------
## 4. Survival subset (drop NA OS and empty groups)
## ------------------------------------------------------------------
if (!all(c("os_time", "os_event") %in% names(dt))) {
  stop("Columns 'os_time' and/or 'os_event' not found in core table.")
}

surv_dt <- dt[!is.na(os_time) & !is.na(os_event) & !is.na(immune_subtype)]
cat("Survival subset dimensions:", paste(dim(surv_dt), collapse = " x "), "\n")

cat("Immune subtype counts in survival subset:\n")
print(table(surv_dt$immune_subtype, useNA = "ifany"))

## Drop subtypes with 0 samples (e.g., LAR)
surv_dt[, immune_subtype := droplevels(immune_subtype)]
subtypes_present <- levels(surv_dt$immune_subtype)
cat("Subtypes present in survival dataset:",
    paste(subtypes_present, collapse = ", "), "\n")

if (length(subtypes_present) < 2) {
  stop("Fewer than 2 immune subtypes present with survival data; cannot run KM.")
}

## ------------------------------------------------------------------
## 5. KM curves by immune_subtype
## ------------------------------------------------------------------
km_out_main <- file.path(
  out_dir,
  "km_immune_subtypes_TNBC_like_TCGA.png"
)

fit <- survfit(Surv(os_time, os_event) ~ immune_subtype, data = surv_dt)

base_palette <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
pal <- base_palette[seq_along(subtypes_present)]

p_km <- ggsurvplot(
  fit,
  data           = surv_dt,
  risk.table     = TRUE,
  pval           = TRUE,
  conf.int       = FALSE,
  legend.title   = "Immune subtype",
  legend.labs    = subtypes_present,
  xlab           = "Time (days)",
  ylab           = "Overall survival probability",
  ggtheme        = theme_minimal(base_size = 13),
  palette        = pal
)

## Save only the KM curve (plot) – avoid adding risk table object
ggsave(km_out_main, p_km$plot, width = 8, height = 6, dpi = 150)
cat("Saved KM plot to:", km_out_main, "\n")

## ------------------------------------------------------------------
## 6. Cox model: immune_subtype (categorical)
## ------------------------------------------------------------------
cox_out <- file.path(
  out_dir,
  "cox_immune_subtypes_TNBC_like_TCGA.tsv"
)

## Use BLIS as reference if present, otherwise first level
if ("BLIS" %in% subtypes_present) {
  surv_dt[, immune_subtype := relevel(immune_subtype, ref = "BLIS")]
} else {
  surv_dt[, immune_subtype := relevel(immune_subtype, ref = subtypes_present[1])]
}

cox_fit <- coxph(Surv(os_time, os_event) ~ immune_subtype, data = surv_dt)
cox_sum <- summary(cox_fit)

hr     <- cox_sum$coef[, "exp(coef)"]
hr_low <- cox_sum$conf.int[, "lower .95"]
hr_high<- cox_sum$conf.int[, "upper .95"]
pvals  <- cox_sum$coef[, "Pr(>|z|)"]

cox_dt <- data.table(
  term   = rownames(cox_sum$coef),
  HR     = as.numeric(hr),
  HR_low = as.numeric(hr_low),
  HR_high= as.numeric(hr_high),
  pvalue = as.numeric(pvals)
)

fwrite(cox_dt, cox_out, sep = "\t")
cat("Saved Cox results to:", cox_out, "\n")

cat("=== DONE Phase II Immune Subtype Summary & Survival ===\n")
