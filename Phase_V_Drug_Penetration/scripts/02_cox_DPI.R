#!/usr/bin/env Rscript

cat("=== Phase V – DPI Survival Models (TNBC-like TCGA) ===\n")

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(survminer)
  library(ggplot2)
})

in_file  <- "results/drug_penetration/DPI_TNBC_like_TCGA.tsv"
out_cox  <- "results/drug_penetration/cox_DPI_TNBC_like_TCGA.tsv"
out_km   <- "results/drug_penetration/KM_DPI_tertiles_TNBC_like_TCGA.png"

if (!file.exists(in_file)) {
  stop("Input file not found: ", in_file)
}

dt <- fread(in_file)
cat("Loaded DPI table: ", nrow(dt), "samples x", ncol(dt), "columns\n")
cat("Columns:\n")
print(names(dt))

## --- Basic safety checks ----------------------------------------------------
needed <- c("os_time", "os_event")
missing_cols <- setdiff(needed, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing required OS columns: ", paste(missing_cols, collapse = ", "))
}

dt[, os_time  := as.numeric(os_time)]
dt[, os_event := as.integer(os_event)]

cat("OS summary:\n")
print(summary(dt$os_time))
cat("Events (os_event==1):", sum(dt$os_event == 1, na.rm = TRUE), "\n")

## --- Surv object ------------------------------------------------------------
surv_obj <- with(dt, Surv(os_time, os_event))

## --- Detect DPI columns -----------------------------------------------------
dpi_cont_candidates  <- intersect(c("DPI_index", "DPI_score", "DPI", "DPI_composite"), names(dt))
dpi_group_candidates <- intersect(c("DPI_tertile", "DPI_group", "DPI_category"), names(dt))

if (length(dpi_cont_candidates) == 0) {
  stop("Could not find a continuous DPI column. Looked for: DPI_index, DPI_score, DPI, DPI_composite")
}
dpi_cont_col <- dpi_cont_candidates[1]

if (length(dpi_group_candidates) == 0) {
  cat("No categorical DPI column found; creating DPI_tertile from ", dpi_cont_col, "\n", sep = "")
  q <- quantile(dt[[dpi_cont_col]], probs = c(1/3, 2/3), na.rm = TRUE)
  dt[, DPI_tertile := cut(
    get(dpi_cont_col),
    breaks = c(-Inf, q[1], q[2], Inf),
    labels = c("Low", "Mid", "High"),
    include.lowest = TRUE
  )]
  dpi_group_col <- "DPI_tertile"
} else {
  dpi_group_col <- dpi_group_candidates[1]
}

cat("Using continuous DPI column: ", dpi_cont_col, "\n", sep = "")
cat("Using categorical DPI column: ", dpi_group_col, "\n", sep = "")

## --- Cox model 1: continuous DPI -------------------------------------------
form1 <- as.formula(paste("Surv(os_time, os_event) ~", dpi_cont_col))
fit1  <- coxph(form1, data = dt)
s1    <- summary(fit1)

m1_hr   <- as.numeric(s1$coef[1, "exp(coef)"])
m1_low  <- as.numeric(s1$conf.int[1, "lower .95"])
m1_high <- as.numeric(s1$conf.int[1, "upper .95"])
m1_p    <- as.numeric(s1$coef[1, "Pr(>|z|)"])

res1 <- data.table(
  model  = "M1_DPI_continuous",
  term   = dpi_cont_col,
  HR     = m1_hr,
  HR_low = m1_low,
  HR_high= m1_high,
  pvalue = m1_p
)

## --- Cox model 2: tertiles vs Mid -------------------------------------------
dt[[dpi_group_col]] <- factor(dt[[dpi_group_col]], levels = c("Mid", "Low", "High"))

form2 <- as.formula(paste("Surv(os_time, os_event) ~", dpi_group_col))
fit2  <- coxph(form2, data = dt)
s2    <- summary(fit2)

coef2 <- s2$coef      # matrix
ci2   <- s2$conf.int  # matrix

terms_raw <- rownames(coef2)

map_term <- function(x) {
  if (grepl("Low",  x)) return(paste0(dpi_group_col, "Low"))
  if (grepl("High", x)) return(paste0(dpi_group_col, "High"))
  return(x)
}

res2 <- data.table(
  model  = "M2_DPI_tertiles_vs_mid",
  term   = vapply(terms_raw, map_term, character(1)),
  HR     = as.numeric(ci2[, "exp(coef)"]),
  HR_low = as.numeric(ci2[, "lower .95"]),
  HR_high= as.numeric(ci2[, "upper .95"]),
  pvalue = as.numeric(coef2[, "Pr(>|z|)"])
)

## --- Save Cox results -------------------------------------------------------
res <- rbind(res1, res2, fill = TRUE)

dir.create(dirname(out_cox), showWarnings = FALSE, recursive = TRUE)
fwrite(res, out_cox, sep = "\t")
cat("Saved Cox results to:", out_cox, "\n")

## --- KM plot: DPI tertiles --------------------------------------------------
dt_km <- dt[!is.na(os_time) & !is.na(os_event) & !is.na(get(dpi_group_col))]

cat("KM subset rows:", nrow(dt_km), "\n")

fit_km <- survfit(as.formula(paste("Surv(os_time, os_event) ~", dpi_group_col)), data = dt_km)

p_km <- ggsurvplot(
  fit_km,
  data            = dt_km,
  risk.table      = TRUE,
  pval            = TRUE,
  legend.title    = dpi_group_col,
  legend.labs     = levels(dt_km[[dpi_group_col]]),
  palette         = c("#1b9e77", "#7570b3", "#d95f02"),
  xlab            = "Days",
  ylab            = "Overall survival probability",
  ggtheme         = theme_minimal()
)

ggsave(out_km, p_km$plot, width = 5, height = 4, dpi = 300)
cat("Saved KM plot to:", out_km, "\n")

cat("=== DONE DPI Cox models ===\n")