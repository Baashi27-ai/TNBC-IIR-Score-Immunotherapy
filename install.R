#!/usr/bin/env Rscript
# =============================================================================
# R Dependencies for TNBC IIR Score Pipeline
# =============================================================================

# Install required packages
packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "survival",
  "survminer",
  "ggplot2",
  "ggpubr",
  "xCell",
  "openxlsx",
  "knitr",
  "rmarkdown"
)

# Install missing packages
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

cat("✅ All R packages installed successfully.\n")
