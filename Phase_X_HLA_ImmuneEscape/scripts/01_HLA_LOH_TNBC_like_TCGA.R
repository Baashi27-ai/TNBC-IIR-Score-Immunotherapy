#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

cat("=== Phase X – HLA-I LOH Detection (TNBC-like TCGA) ===\n")

## Always work from project root
setwd("/mnt/c/TNBC_project")

## 1) Locate TNBC MAF file (same as earlier TMB-from-MAF step)
maf_candidates <- c(
  "Biomarker_Verification/inputs/genomic/TNBC_only.maf",
  "data_proc/TNBC_only.maf"
)

maf_file <- NULL
for (f in maf_candidates) {
  if (file.exists(f)) {
    maf_file <- f
    break
  }
}

if (is.null(maf_file)) {
  stop("Could not find TNBC_only.maf in expected locations. Checked:\n",
       paste(maf_candidates, collapse = "\n"))
}

cat("Using MAF file:\n  ", maf_file, "\n")

## 2) Read MAF
maf <- fread(maf_file)
cat("Loaded MAF rows:", nrow(maf), "\n")

## 3) Make a unified submitter_id column
id_col <- NULL
for (col in c("submitter_id", "Tumor_Sample_Barcode", "tumor_sample_barcode",
              "Sample_ID", "sample", "Tumor_Sample")) {
  if (col %in% names(maf)) {
    id_col <- col
    break
  }
}

if (is.null(id_col)) {
  stop("Could not find any standard sample ID column in MAF.\n",
       "Checked: submitter_id, Tumor_Sample_Barcode, Sample_ID, sample, Tumor_Sample")
}

cat("Using sample column for ID:", id_col, "\n")

if (id_col == "submitter_id") {
  maf[, submitter_id := get(id_col)]
} else {
  ## Derive 12-char TCGA submitter ID: TCGA-XX-XXXX
  maf[, submitter_id := substr(get(id_col), 1L, 12L)]
}

## 4) Filter for class-I HLA genes
hla_genes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G")

if (!"Hugo_Symbol" %in% names(maf)) {
  stop("MAF does not contain 'Hugo_Symbol' column; cannot select HLA genes.")
}

maf_hla <- maf[Hugo_Symbol %in% hla_genes]

cat("HLA variants retained:", nrow(maf_hla), "\n")

out_dir <- "Immune_Spatial_Immunotherapy/Phase_X_HLA_LOH/results/hla_loh"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, "HLA_LOH_proxy_TNBC_like_TCGA.tsv")

if (nrow(maf_hla) == 0L) {
  warning("No HLA variants found in MAF; writing empty result table.")
  fwrite(data.table(), file = out_file, sep = "\t")
  cat("Saved EMPTY HLA LOH proxy table to:\n  ", out_file, "\n")
  quit(save = "no")
}

## 5) Mark “damaging” variants (LOH-like hits)
damaging_classes <- c(
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "Nonsense_Mutation",
  "Splice_Site",
  "Translation_Start_Site",
  "Nonstop_Mutation"
)

if (!"Variant_Classification" %in% names(maf_hla)) {
  stop("MAF does not contain 'Variant_Classification'; cannot classify damaging variants.")
}

maf_hla[, damaging := fifelse(Variant_Classification %in% damaging_classes, 1L, 0L)]

## 6) Aggregate per patient + HLA gene
dt_gene <- maf_hla[
  ,
  .(
    n_variants = .N,
    n_damaging = sum(damaging, na.rm = TRUE)
  ),
  by = .(submitter_id, Hugo_Symbol)
]

## 7) Wide gene-level table: columns per HLA gene
dt_wide <- dcast(
  dt_gene,
  submitter_id ~ Hugo_Symbol,
  value.var = "n_damaging",
  fill = 0
)

## 8) Patient-level summary metrics
dt_summary <- dt_gene[
  ,
  .(
    HLA_total_variants    = sum(n_variants),
    HLA_damaging_variants = sum(n_damaging),
    HLA_damaging_prop     = ifelse(sum(n_variants) > 0,
                                   sum(n_damaging) / sum(n_variants),
                                   NA_real_)
  ),
  by = submitter_id
]

## 9) Merge summary + per-gene
hla_res <- merge(
  dt_summary,
  dt_wide,
  by = "submitter_id",
  all.x = TRUE
)

## 10) Define an LOH-proxy group by tertiles of damaging hits
valid <- hla_res[!is.na(HLA_damaging_variants)]

if (nrow(valid) > 0L) {
  qs <- quantile(valid$HLA_damaging_variants,
                 probs = c(1/3, 2/3),
                 na.rm = TRUE)

  hla_res[
    ,
    HLA_LOH_proxy_group := fifelse(
      is.na(HLA_damaging_variants), NA_character_,
      fifelse(
        HLA_damaging_variants <= qs[1], "Low",
        fifelse(HLA_damaging_variants <= qs[2], "Mid", "High")
      )
    )
  ]
} else {
  hla_res[, HLA_LOH_proxy_group := NA_character_]
}

## 11) Save
fwrite(hla_res, out_file, sep = "\t")
cat("Saved HLA LOH proxy table to:\n  ", out_file, "\n")

cat("Preview:\n")
print(head(hla_res))

cat("=== DONE Phase X HLA-I LOH proxy ===\n")