#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(survminer)
})

message("=== Phase I – Immune Survival Modeling (TNBC-like TCGA) ===")

infile <- "results/immune_scores/immune_hot_cold_TNBC_like_TCGA_OS.tsv"
out_km  <- "results/immune_scores/km_immune_hot_vs_cold.png"
out_km3 <- "results/immune_scores/km_immune_3groups.png"
out_cox <- "results/immune_scores/cox_immune_results.tsv"

# --------------------------
# 1) Load data
# --------------------------
dt <- fread(infile)
message("Loaded: ", nrow(dt), " samples")

# Clean OS
dt[, os_time := as.numeric(os_time)]
dt[, os_event := as.integer(os_event)]

dt <- dt[!is.na(os_time) & !is.na(os_event)]

# --------------------------
# 2) KM: Hot vs Cold
# --------------------------
dt2 <- dt[immune_group %in% c("Hot", "Cold")]

fit2 <- survfit(Surv(os_time, os_event) ~ immune_group, data = dt2)

p1 <- ggsurvplot(
  fit2,
  data = dt2,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  legend.title = "Immune Group",
  legend.labs = c("Cold", "Hot"),
  palette = c("#1f77b4", "#d62728")
)

ggsave(out_km, p1$plot, width = 6, height = 5)
message("Saved: ", out_km)

# --------------------------
# 3) KM: 3 groups
# --------------------------
fit3 <- survfit(Surv(os_time, os_event) ~ immune_group, data = dt)

p2 <- ggsurvplot(
  fit3,
  data = dt,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  legend.title = "Immune Group",
  palette = c("Cold" = "#1f77b4", "Intermediate" = "#2ca02c", "Hot" = "#d62728")
)

ggsave(out_km3, p2$plot, width = 6, height = 5)
message("Saved: ", out_km3)

# --------------------------
# 4) Cox models
# --------------------------

# continuous ImmuneScore
cox1 <- coxph(Surv(os_time, os_event) ~ ImmuneScore, data = dt)
sum1 <- summary(cox1)

# categorical immune_group
cox2 <- coxph(Surv(os_time, os_event) ~ immune_group, data = dt)
sum2 <- summary(cox2)

# Save results
res <- data.table(
  model = c("ImmuneScore_continuous", "immune_group_categorical"),
  HR = c(sum1$coef[1,2], NA),
  pvalue = c(sum1$coef[1,5], sum2$wald["pvalue"]),
  HR_low = c(sum1$conf.int[,"lower .95"], NA),
  HR_high = c(sum1$conf.int[,"upper .95"], NA)
)

fwrite(res, out_cox, sep = "\t")
message("Saved COX results to: ", out_cox)

message("=== DONE Immune KM & COX ===")
