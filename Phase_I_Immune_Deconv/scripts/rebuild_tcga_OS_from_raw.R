#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

message("=== Rebuild real TCGA OS file from raw clinical ===")

# 1) Try richer clinical file first, then fallback
path1 <- "data_raw/tcga/clinical/tcga_clinical_with_subtype_proxy.csv"
path2 <- "data_raw/tcga/clinical/tcga_clinical.csv"

if (file.exists(path1)) {
  clin_raw_path <- path1
} else if (file.exists(path2)) {
  clin_raw_path <- path2
} else {
  stop("Could not find tcga_clinical_with_subtype_proxy.csv or tcga_clinical.csv under data_raw/tcga/clinical/")
}

message("Using raw clinical file: ", clin_raw_path)

clin <- fread(clin_raw_path)
message("Raw clinical dimensions: ", paste(dim(clin), collapse = " x "))
message("Clinical columns: ", paste(names(clin), collapse = ", "))

## 2) Choose an ID column
id_candidates <- c(
  "bcr_patient_barcode",
  "submitter_id",
  "case_submitter_id",
  "PATIENT_ID",
  "patient_id"
)

id_col <- id_candidates[id_candidates %in% names(clin)][1]

if (is.na(id_col)) {
  stop("Could not find any suitable ID column (bcr_patient_barcode / submitter_id / case_submitter_id / PATIENT_ID / patient_id)")
}

message("Using ID column: ", id_col)

# Standardize: submitter_id (case-level, 12 chars)
clin[, submitter_id := substr(get(id_col), 1, 12)]
clin[, case_id := submitter_id]

## 3) Build OS time and event

# First, try to use direct OS columns (your file has OS_time_days and OS_status)
time_candidates <- c("OS_time_days", "OS.time", "os_time", "OS_days", "os_time_days")
event_candidates <- c("OS_status", "OS", "os_event", "OS_EVENT")

time_col <- time_candidates[time_candidates %in% names(clin)][1]
event_col <- event_candidates[event_candidates %in% names(clin)][1]

have_direct_OS <- !is.na(time_col) && !is.na(event_col)

if (have_direct_OS) {
  message("Found direct OS columns: ", time_col, " / ", event_col)
  clin[, os_time := as.numeric(get(time_col))]

  # OS_status is usually 1/0 or 'DECEASED'/'ALIVE'
  raw_event <- clin[[event_col]]

  if (is.character(raw_event)) {
    vs <- toupper(trimws(raw_event))
    clin[, os_event := fifelse(vs %in% c("1", "DECEASED", "DEAD"), 1L,
                         fifelse(vs %in% c("0", "ALIVE"), 0L, NA_integer_))]
  } else {
    clin[, os_event := as.integer(raw_event)]
  }

} else {
  message("No direct OS columns found. Attempting to build from days_to_death / follow_up")

  d_death_cols <- c("days_to_death", "DAYS_TO_DEATH")
  d_fu_cols    <- c("days_to_last_follow_up", "DAYS_TO_LAST_FOLLOWUP")
  vital_cols   <- c("vital_status", "VITAL_STATUS")

  d_death_col <- d_death_cols[d_death_cols %in% names(clin)][1]
  d_fu_col    <- d_fu_cols[d_fu_cols %in% names(clin)][1]
  vital_col   <- vital_cols[vital_cols %in% names(clin)][1]

  if (is.na(d_fu_col) && is.na(d_death_col)) {
    stop("No days_to_death / days_to_last_follow_up columns found; cannot construct OS.")
  }

  if (!is.na(d_death_col)) {
    clin[, d_death := suppressWarnings(as.numeric(get(d_death_col)))]
  } else {
    clin[, d_death := NA_real_]
  }

  if (!is.na(d_fu_col)) {
    clin[, d_fu := suppressWarnings(as.numeric(get(d_fu_col)))]
  } else {
    clin[, d_fu := NA_real_]
  }

  if (!is.na(vital_col)) {
    vs <- toupper(trimws(clin[[vital_col]]))
    clin[, os_event := fifelse(vs %in% c("DECEASED", "DEAD"), 1L,
                         fifelse(vs == "ALIVE", 0L, NA_integer_))]
    clin[is.na(os_event) & !is.na(d_death) & d_death > 0, os_event := 1L]
    clin[is.na(os_event) & is.na(d_death), os_event := 0L]
  } else {
    clin[, os_event := fifelse(!is.na(d_death) & d_death > 0, 1L, 0L)]
  }

  clin[, os_time := fifelse(!is.na(d_death), d_death, d_fu)]
}

## 4) Keep a compact OS table
keep_cols <- c("case_id", "submitter_id", "os_time", "os_event")
keep_cols <- keep_cols[keep_cols %in% names(clin)]

clin_os <- unique(clin[, ..keep_cols])
setorder(clin_os, submitter_id)

message("Final OS table dimensions: ", paste(dim(clin_os), collapse = " x "))
message("Non-NA os_time rows: ", sum(!is.na(clin_os$os_time)))
message("Non-NA os_event rows:", sum(!is.na(clin_os$os_event)))

## 5) Write to both locations
out1 <- "Biomarker_Verification/inputs/external_cohorts/tcga_brca_clinical_case12_OSbuild.csv"
out2 <- "Immune_Spatial_Immunotherapy/Phase_I_Immune_Deconv/inputs/clinical/tcga_brca_clinical_case12_OSbuild.csv"

dir.create(dirname(out1), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out2), recursive = TRUE, showWarnings = FALSE)

fwrite(clin_os, out1)
fwrite(clin_os, out2)

message("Wrote OS table to: ", out1)
message("Wrote OS table to: ", out2)

message("Preview:")
print(head(clin_os, 10))
message("=== DONE rebuilding OS file ===")
