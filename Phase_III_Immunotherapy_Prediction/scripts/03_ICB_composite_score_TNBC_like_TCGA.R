#!/usr/bin/env Rscript

cat("=== Phase III – ICB Composite Cox Models (TNBC-like TCGA) ===\n")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(survival)
})

# ---------- Input / output paths ----------
pd1_file  <- "results/immunotherapy/PD1_PDL1_signature_TNBC_like_TCGA.tsv"
tmb_file  <- "results/immunotherapy/TMB_TNBC_like_TCGA.tsv"
out_dir   <- "results/immunotherapy"
out_cox   <- file.path(out_dir, "Cox_ICB_models_TNBC_like_TCGA.tsv")

# ---------- Load PD1/PDL1 signature + OS ----------
cat("Reading PD1/PDL1 signature file from:\n  ", pd1_file, "\n", sep = "")
pd1_dt <- fread(pd1_file)

required_cols <- c("submitter_id", "os_time", "os_event", "PD1_PDL1_signature")
missing_pd1   <- setdiff(required_cols, names(pd1_dt))
if (length(missing_pd1) > 0) {
  stop("PD1/PDL1 file missing required columns: ", paste(missing_pd1, collapse = ", "))
}

cat("PD1 table dimensions:", paste(dim(pd1_dt), collapse = " x "), "\n")
cat("PD1 summary (signature):\n")
print(summary(pd1_dt$PD1_PDL1_signature))

# ---------- Load TMB (may be placeholder with NAs) ----------
cat("Reading TMB file from:\n  ", tmb_file, "\n", sep = "")
tmb_dt <- fread(tmb_file)

if (!"submitter_id" %in% names(tmb_dt)) {
  stop("TMB file does not contain column 'submitter_id'")
}

if (!"TMB" %in% names(tmb_dt)) {
  cat("WARNING: TMB column not found in TMB file. Creating TMB = NA.\n")
  tmb_dt[, TMB := NA_real_]
}

cat("TMB table dimensions:", paste(dim(tmb_dt), collapse = " x "), "\n")
cat("Non-NA TMB values:", sum(!is.na(tmb_dt$TMB)), "\n")

# ---------- Merge PD1 + TMB ----------
full_dt <- pd1_dt %>%
  left_join(
    tmb_dt %>% select(submitter_id, TMB),
    by = "submitter_id"
  )

cat("Merged ICB dataset dimensions:", paste(dim(full_dt), collapse = " x "), "\n")

# Ensure os_event is numeric (0/1)
full_dt <- full_dt %>%
  mutate(
    os_event = as.integer(os_event)
  )

# ---------- Define PD1/PDL1 high vs low (median split) ----------
med_pd1 <- median(full_dt$PD1_PDL1_signature, na.rm = TRUE)
cat("Median PD1/PDL1 signature:", med_pd1, "\n")

full_dt <- full_dt %>%
  mutate(
    PD1_group = ifelse(PD1_PDL1_signature >= med_pd1, "High", "Low")
  )

cat("PD1 high/low group counts:\n")
print(table(full_dt$PD1_group))

# ---------- Helper: extract Cox results ----------
extract_cox <- function(fit, model_name) {
  s <- summary(fit)
  coefs <- as.data.frame(s$coef)
  cis   <- as.data.frame(s$conf.int)

  # Column names in conf.int are like: "exp(coef)", "exp(-coef)", "lower .95", "upper .95"
  hr      <- cis[, "exp(coef)"]
  hr_low  <- cis[, "lower .95"]
  hr_high <- cis[, "upper .95"]
  pval    <- coefs[, "Pr(>|z|)"]
  terms   <- rownames(coefs)

  data.frame(
    model   = model_name,
    term    = terms,
    HR      = hr,
    HR_low  = hr_low,
    HR_high = hr_high,
    pvalue  = pval,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

cox_rows <- list()

# ---------- Model 1: PD1/PDL1 signature (continuous) ----------
cat("\n--- Model 1: PD1/PDL1 continuous ---\n")
d1 <- full_dt %>%
  filter(!is.na(os_time), !is.na(os_event), !is.na(PD1_PDL1_signature))

cat("Model 1 N:", nrow(d1), "events:", sum(d1$os_event, na.rm = TRUE), "\n")

if (nrow(d1) > 10 && sum(d1$os_event, na.rm = TRUE) >= 5) {
  fit1 <- coxph(Surv(os_time, os_event) ~ PD1_PDL1_signature, data = d1)
  print(summary(fit1))
  c1 <- extract_cox(fit1, "M1_PD1_continuous")
  cox_rows[["M1"]] <- c1
} else {
  cat("Skipping Model 1: not enough samples/events.\n")
}

# ---------- Model 2: PD1 group (High vs Low) ----------
cat("\n--- Model 2: PD1 high vs low ---\n")
d2 <- full_dt %>%
  filter(!is.na(os_time), !is.na(os_event), !is.na(PD1_group))

cat("Model 2 N:", nrow(d2), "events:", sum(d2$os_event, na.rm = TRUE), "\n")
cat("Group counts:\n")
print(table(d2$PD1_group))

if (nrow(d2) > 10 && sum(d2$os_event, na.rm = TRUE) >= 5) {
  d2$PD1_group <- factor(d2$PD1_group, levels = c("Low", "High"))
  fit2 <- coxph(Surv(os_time, os_event) ~ PD1_group, data = d2)
  print(summary(fit2))
  c2 <- extract_cox(fit2, "M2_PD1_high_vs_low")
  cox_rows[["M2"]] <- c2
} else {
  cat("Skipping Model 2: not enough samples/events.\n")
}

# ---------- Model 3: TMB continuous (if available) ----------
cat("\n--- Model 3: TMB continuous ---\n")
d3 <- full_dt %>%
  filter(!is.na(os_time), !is.na(os_event), !is.na(TMB))

cat("Model 3 N:", nrow(d3), "events:", sum(d3$os_event, na.rm = TRUE), "\n")

if (nrow(d3) > 10 && sum(d3$os_event, na.rm = TRUE) >= 5) {
  fit3 <- coxph(Surv(os_time, os_event) ~ TMB, data = d3)
  print(summary(fit3))
  c3 <- extract_cox(fit3, "M3_TMB_continuous")
  cox_rows[["M3"]] <- c3
} else {
  cat("Skipping Model 3: not enough non-NA TMB / events.\n")
}

# ---------- Model 4: PD1 + TMB (joint model, if TMB exists) ----------
cat("\n--- Model 4: PD1 + TMB joint ---\n")
d4 <- full_dt %>%
  filter(!is.na(os_time), !is.na(os_event),
         !is.na(PD1_PDL1_signature), !is.na(TMB))

cat("Model 4 N:", nrow(d4), "events:", sum(d4$os_event, na.rm = TRUE), "\n")

if (nrow(d4) > 10 && sum(d4$os_event, na.rm = TRUE) >= 5) {
  fit4 <- coxph(Surv(os_time, os_event) ~ PD1_PDL1_signature + TMB, data = d4)
  print(summary(fit4))
  c4 <- extract_cox(fit4, "M4_PD1_plus_TMB")
  cox_rows[["M4"]] <- c4
} else {
  cat("Skipping Model 4: not enough non-NA TMB / events.\n")
}

# ---------- Combine & write ----------
if (length(cox_rows) == 0) {
  cat("\nWARNING: No Cox models were successfully fitted. No output file will be created.\n")
} else {
  cox_all <- bind_rows(cox_rows)
  fwrite(as.data.table(cox_all), out_cox, sep = "\t")
  cat("\nWrote Cox ICB models to:\n  ", out_cox, "\n", sep = "")
}

cat("=== DONE Phase III ICB composite Cox ===\n")
