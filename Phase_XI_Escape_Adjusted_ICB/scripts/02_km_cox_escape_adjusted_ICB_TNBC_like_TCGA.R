#!/usr/bin/env Rscript

cat("=== Phase XI – KM & Cox for Escape-Adjusted ICB (TNBC-like TCGA) ===\n")

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(survminer)
  library(ggplot2)
})

in_file  <- "results/escape_adjusted/escape_adjusted_ICB_TNBC_like_TCGA.tsv"
out_dir  <- "results/escape_adjusted"

if (!file.exists(in_file)) {
  stop("Input escape-adjusted table not found at: ", in_file)
}

dt <- fread(in_file)
cat("Loaded escape-adjusted table:\n  ", in_file, "\n")
cat("Dimensions:", nrow(dt), "x", ncol(dt), "\n")

## 1) Survival subset
dt[, os_time  := as.numeric(os_time)]
dt[, os_event := as.numeric(os_event)]

dt_surv <- dt[!is.na(os_time) & !is.na(os_event) & os_time > 0]

cat("Survival subset rows:", nrow(dt_surv), "\n")
cat("OS summary:\n")
print(summary(dt_surv$os_time))

cat("Events (os_event == 1):", dt_surv[os_event == 1, .N], "\n")

# Ensure group is factor with sensible order
if (!"escape_adjusted_ICB_group" %in% colnames(dt_surv)) {
  stop("Column 'escape_adjusted_ICB_group' not found.")
}

dt_surv[, escape_adjusted_ICB_group := factor(
  escape_adjusted_ICB_group,
  levels = c("Escape_adjusted_Poor_ICB_candidate",
             "Escape_adjusted_Intermediate_candidate",
             "Escape_adjusted_High_ICB_candidate")
)]

cat("Group counts (survival subset):\n")
print(dt_surv[, .N, by = escape_adjusted_ICB_group])

## Build Surv object
surv_obj <- with(dt_surv, Surv(os_time, os_event))

## 2) Cox model – continuous escape-adjusted score
cox_cont <- coxph(surv_obj ~ escape_adjusted_ICB_score, data = dt_surv)
summary_cont <- summary(cox_cont)

HR_cont      <- summary_cont$coef[1, "exp(coef)"]
HR_low_cont  <- summary_cont$conf.int[1, "lower .95"]
HR_high_cont <- summary_cont$conf.int[1, "upper .95"]
p_cont       <- summary_cont$coef[1, "Pr(>|z|)"]

## 3) Cox model – categorical groups (ref = Intermediate)
dt_surv[, escape_group_relevel := relevel(
  escape_adjusted_ICB_group,
  ref = "Escape_adjusted_Intermediate_candidate"
)]

cox_cat <- coxph(surv_obj ~ escape_group_relevel, data = dt_surv)
summary_cat <- summary(cox_cat)

cox_rows <- list()

# Continuous row
cox_rows[[length(cox_rows) + 1]] <- data.table(
  model  = "M1_escape_adjusted_continuous",
  term   = "escape_adjusted_ICB_score",
  HR     = HR_cont,
  HR_low = HR_low_cont,
  HR_high= HR_high_cont,
  pvalue = p_cont
)

# Categorical rows
coef_cat <- summary_cat$coef
ci_cat   <- summary_cat$conf.int

for (i in seq_len(nrow(coef_cat))) {
  term_name <- rownames(coef_cat)[i]
  HR        <- ci_cat[i, "exp(coef)"]
  HR_low    <- ci_cat[i, "lower .95"]
  HR_high   <- ci_cat[i, "upper .95"]
  pval      <- coef_cat[i, "Pr(>|z|)"]

  cox_rows[[length(cox_rows) + 1]] <- data.table(
    model  = "M2_escape_groups_vs_Intermediate",
    term   = term_name,
    HR     = HR,
    HR_low = HR_low,
    HR_high= HR_high,
    pvalue = pval
  )
}

cox_dt <- rbindlist(cox_rows, use.names = TRUE, fill = TRUE)

# Save Cox results
cox_file <- file.path(out_dir, "cox_escape_adjusted_ICB_TNBC_like_TCGA.tsv")
fwrite(cox_dt, cox_file, sep = "\t")
cat("Saved Cox results to:\n  ", cox_file, "\n")
print(cox_dt)

## 4) KM curves for escape-adjusted groups
fit_km <- survfit(surv_obj ~ escape_adjusted_ICB_group, data = dt_surv)

km_plot <- ggsurvplot(
  fit_km,
  data        = dt_surv,
  risk.table  = TRUE,
  pval        = TRUE,
  conf.int    = FALSE,
  legend.title= "Escape-adjusted ICB group",
  legend.labs = c("Poor", "Intermediate", "High"),
  xlab        = "Overall survival (days)",
  ylab        = "Survival probability"
)

# Save main KM plot
km_file <- file.path(out_dir, "KM_escape_adjusted_ICB_groups_TNBC_like_TCGA.png")
ggsave(filename = km_file, plot = km_plot$plot, width = 7, height = 6, dpi = 300)
cat("Saved KM plot to:\n  ", km_file, "\n")

# Save KM + risk table as separate file (optional)
km_rt_file <- file.path(out_dir, "KM_escape_adjusted_ICB_groups_TNBC_like_TCGA_risktable.png")
ggsave(filename = km_rt_file, plot = km_plot$table, width = 7, height = 4, dpi = 300)
cat("Saved KM risk-table plot to:\n  ", km_rt_file, "\n")

cat("=== DONE Phase XI – KM & Cox for Escape-Adjusted ICB ===\n")