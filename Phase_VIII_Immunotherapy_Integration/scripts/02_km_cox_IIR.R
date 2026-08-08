#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(survminer)
  library(ggplot2)
})

cat("=== Phase VIII – KM & Cox for Integrated Immunotherapy Readiness (TNBC-like TCGA) ===\n")

## ----------------------------
## 1) Paths
## ----------------------------

in_file  <- "results/integration/IIR_table_TNBC_like_TCGA.tsv"
out_dir  <- "results/integration"
out_cox  <- file.path(out_dir, "cox_IIR_TNBC_like_TCGA.tsv")
out_km_p <- file.path(out_dir, "KM_IIR_groups_TNBC_like_TCGA.png")
out_km_r <- file.path(out_dir, "KM_IIR_groups_TNBC_like_TCGA_risktable.png")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

## ----------------------------
## 2) Load integrated table
## ----------------------------

cat("Reading integrated table from:\n   ", in_file, "\n")
dt <- fread(in_file)
cat("Original rows:", nrow(dt), "\n")
cat("Columns:\n")
print(colnames(dt))

## ----------------------------
## 3) Check key columns
## ----------------------------

required_cols <- c("submitter_id", "os_time", "os_event", "IIR_score", "IIR_group")
missing <- setdiff(required_cols, colnames(dt))
if (length(missing) > 0) {
  stop("Missing required columns in IIR table: ", paste(missing, collapse = ", "))
}

dt[, IIR_group := factor(IIR_group,
                         levels = c("Poor_ICB_ready", "Intermediate", "High_ICB_ready"))]

cat("IIR_group counts (all rows):\n")
print(table(dt$IIR_group, useNA = "ifany"))

## ----------------------------
## 4) Subset for survival analysis
## ----------------------------

dt_surv <- dt[!is.na(os_time) & !is.na(os_event) & !is.na(IIR_group)]

cat("Survival subset rows:", nrow(dt_surv), "\n")
cat("OS summary in survival subset:\n")
print(summary(dt_surv$os_time))
cat("Events (os_event == 1):", sum(dt_surv$os_event == 1, na.rm = TRUE), "\n")
cat("IIR_group counts (survival subset):\n")
print(table(dt_surv$IIR_group, useNA = "ifany"))

if (nrow(dt_surv) < 20) {
  warning("Very small survival subset – interpretation will be unstable.")
}

## ----------------------------
## 5) Safe Cox extractor
## ----------------------------

extract_cox <- function(fit, model_name) {
  s    <- summary(fit)
  coefs <- as.data.frame(s$coefficients)
  cis   <- as.data.frame(s$conf.int)

  # Find columns by pattern instead of using weird names like "lower .95"
  hr_col   <- grep("exp\\(coef\\)", colnames(cis), value = FALSE)
  low_col  <- grep("lower",        colnames(cis), value = FALSE)
  high_col <- grep("upper",        colnames(cis), value = FALSE)
  p_col    <- grep("^Pr",          colnames(coefs), value = FALSE)

  res <- data.table(
    model   = model_name,
    term    = rownames(coefs),
    HR      = unname(cis[, hr_col]),
    HR_low  = unname(cis[, low_col]),
    HR_high = unname(cis[, high_col]),
    pvalue  = unname(coefs[, p_col])
  )
  res
}

cox_results <- list()

## Model 1: continuous IIR_score
cat("Fitting Cox model: IIR_score (continuous)...\n")
fit1 <- coxph(Surv(os_time, os_event) ~ IIR_score, data = dt_surv)
cox_results[[length(cox_results) + 1]] <- extract_cox(fit1, "M1_IIR_continuous")

## Model 2: categorical IIR_group (Intermediate as reference)
cat("Fitting Cox model: IIR_group (categorical, ref = Intermediate)...\n")
dt_surv[, IIR_group_ref := relevel(IIR_group, ref = "Intermediate")]
fit2 <- coxph(Surv(os_time, os_event) ~ IIR_group_ref, data = dt_surv)
cox_results[[length(cox_results) + 1]] <- extract_cox(fit2, "M2_IIR_groups_vs_Intermediate")

## Optional Model 3: add TMB_tertile if present
if ("TMB_tertile" %in% colnames(dt_surv)) {
  cat("Fitting Cox model: IIR_score + TMB_tertile (optional)...\n")
  dt_surv[, TMB_tertile := factor(TMB_tertile, levels = c("Low", "Mid", "High"))]
  fit3 <- coxph(Surv(os_time, os_event) ~ IIR_score + TMB_tertile, data = dt_surv)
  cox_results[[length(cox_results) + 1]] <- extract_cox(fit3, "M3_IIR_plus_TMBtertiles")
}

cox_dt <- rbindlist(cox_results, use.names = TRUE, fill = TRUE)

cat("Saving Cox results to:\n   ", out_cox, "\n")
fwrite(cox_dt, out_cox, sep = "\t")

## ----------------------------
## 6) KM plot for IIR_group
## ----------------------------

cat("Fitting KM curves by IIR_group...\n")
fit_km <- survfit(Surv(os_time, os_event) ~ IIR_group, data = dt_surv)

cat("Creating KM plot...\n")
p <- ggsurvplot(
  fit_km,
  data         = dt_surv,
  risk.table   = TRUE,
  pval         = TRUE,
  conf.int     = FALSE,
  xlab         = "Time (days)",
  ylab         = "Overall survival probability",
  legend.title = "IIR group",
  legend.labs  = c("Poor ICB-ready", "Intermediate", "High ICB-ready")
)

cat("Saving KM plot to:\n   ", out_km_p, "\n")
ggsave(out_km_p, p$plot, width = 6, height = 5, dpi = 300)

cat("Saving KM risk table to:\n   ", out_km_r, "\n")
ggsave(out_km_r, p$table, width = 6, height = 2.5, dpi = 300)

cat("=== DONE KM & Cox for IIR ===\n")