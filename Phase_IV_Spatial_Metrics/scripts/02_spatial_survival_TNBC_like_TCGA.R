#!/usr/bin/env Rscript

cat("=== Phase IV – Spatial Metrics Survival Modeling (TNBC-like TCGA) ===\n")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(survival)
  library(survminer)
  library(ggplot2)
})

in_file  <- file.path("results", "spatial_metrics",
                      "spatial_metrics_with_clin_TNBC_like_TCGA.tsv")
out_dir  <- file.path("results", "spatial_metrics")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Reading merged spatial + clinical table from:\n  ", in_file, "\n")
dt <- fread(in_file)

cat("Merged table dimensions:", nrow(dt), "x", ncol(dt), "\n")
cat("Columns:", paste(colnames(dt), collapse = ", "), "\n")

# Basic cleaning for survival
surv_dt <- dt %>%
  filter(!is.na(os_time), !is.na(os_event)) %>%
  mutate(
    os_time  = as.numeric(os_time),
    os_event = as.integer(os_event)
  )

cat("Survival subset dimensions:", nrow(surv_dt), "x", ncol(surv_dt), "\n")

# Ensure factor levels for tertiles
surv_dt$Immune_Stroma_tertile <- factor(
  surv_dt$Immune_Stroma_tertile,
  levels = c("Low", "Mid", "High")
)

cat("Immune_Stroma_tertile counts (survival subset):\n")
print(table(surv_dt$Immune_Stroma_tertile, useNA = "ifany"))

# -----------------------------
# Cox models
# -----------------------------
cox_res <- list()

# Model 1 – continuous Immune:Stroma ratio
cox1 <- coxph(Surv(os_time, os_event) ~ Immune_Stroma_ratio, data = surv_dt)
summary1 <- summary(cox1)

cox_res[["M1_ImmuneStroma_continuous"]] <- data.table(
  model  = "M1_ImmuneStroma_continuous",
  term   = "Immune_Stroma_ratio",
  HR     = summary1$coefficients[,"exp(coef)"],
  HR_low = summary1$conf.int[,"lower .95"],
  HR_high= summary1$conf.int[,"upper .95"],
  pvalue = summary1$coefficients[,"Pr(>|z|)"]
)

# Model 2 – tertiles (reference = Mid)
surv_dt$Immune_Stroma_tertile <- relevel(surv_dt$Immune_Stroma_tertile, ref = "Mid")

cox2 <- coxph(Surv(os_time, os_event) ~ Immune_Stroma_tertile, data = surv_dt)
summary2 <- summary(cox2)

for (i in seq_len(nrow(summary2$coefficients))) {
  nm <- rownames(summary2$coefficients)[i]
  cox_res[[paste0("M2_", nm)]] <- data.table(
    model  = "M2_ImmuneStroma_tertiles_vs_mid",
    term   = nm,
    HR     = summary2$coefficients[i, "exp(coef)"],
    HR_low = summary2$conf.int[i, "lower .95"],
    HR_high= summary2$conf.int[i, "upper .95"],
    pvalue = summary2$coefficients[i, "Pr(>|z|)"]
  )
}

# Bind results
cox_out <- rbindlist(cox_res, use.names = TRUE, fill = TRUE)

cox_file <- file.path(out_dir, "cox_spatial_metrics_TNBC_like_TCGA.tsv")
fwrite(cox_out, cox_file, sep = "\t")

cat("Saved Cox results to:\n  ", cox_file, "\n")
cat("Cox results preview:\n")
print(cox_out)

# -----------------------------
# KM plot – High vs Low Immune:Stroma (drop Mid)
# -----------------------------
km_dt <- surv_dt %>%
  filter(Immune_Stroma_tertile %in% c("Low", "High")) %>%
  droplevels()

cat("KM High vs Low subset counts:\n")
print(table(km_dt$Immune_Stroma_tertile, useNA = "ifany"))

fit_km <- survfit(Surv(os_time, os_event) ~ Immune_Stroma_tertile, data = km_dt)

p_km <- ggsurvplot(
  fit_km,
  data        = km_dt,
  risk.table  = TRUE,
  pval        = TRUE,
  conf.int    = FALSE,
  legend.title= "Immune:Stroma",
  legend.labs = c("Low", "High"),
  xlab        = "Time (days)",
  ylab        = "Overall survival probability",
  ggtheme     = theme_bw()
)

km_file <- file.path(out_dir, "KM_ImmuneStroma_High_vs_Low_TNBC_like_TCGA.png")
km_rt   <- file.path(out_dir, "KM_ImmuneStroma_High_vs_Low_TNBC_like_TCGA_risktable.png")

ggsave(km_file, p_km$plot, width = 6, height = 5, dpi = 300)
ggsave(km_rt,   p_km$table + theme_bw(), width = 6, height = 3, dpi = 300)

cat("Saved KM plot to:\n  ", km_file, "\n")
cat("Saved KM risk table to:\n  ", km_rt, "\n")

# -----------------------------
# Optional – Heatmap-style barplot: Immune subtype vs spatial tertile
# -----------------------------
p_ct <- ggplot(surv_dt, aes(x = immune_subtype, fill = Immune_Stroma_tertile)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw() +
  xlab("Immune subtype") +
  ylab("Proportion") +
  ggtitle("Distribution of Immune:Stroma tertiles across immune subtypes")

ct_file <- file.path(out_dir, "stackedbar_ImmuneStroma_by_immune_subtype_TNBC_like_TCGA.png")
ggsave(ct_file, p_ct, width = 6, height = 4, dpi = 300)

cat("Saved subtype x spatial stacked barplot to:\n  ", ct_file, "\n")
cat("=== DONE Phase IV – Spatial survival modeling ===\n")