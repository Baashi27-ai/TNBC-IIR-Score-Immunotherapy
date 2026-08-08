#!/usr/bin/env Rscript

cat("=== Phase VII – SBS96 spectrum + Ageing/APOBEC scores (TNBC-like TCGA) ===\n")

suppressPackageStartupMessages({
  library(data.table)
  library(BiocManager)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings)
})

#---------- Paths ----------
maf_file <- "inputs/maf/TNBC_only.maf"
out_mat  <- "results/signatures/SBS96_mutation_matrix_TNBC_like_TCGA.tsv"
out_sc   <- "results/signatures/signature_scores_TNBC_like_TCGA.tsv"

# Ensure results dir exists
dir.create("results/signatures", showWarnings = FALSE, recursive = TRUE)

cat("Reading MAF from:\n  ", maf_file, "\n")
maf <- fread(maf_file)
cat("MAF rows:", nrow(maf), "\n")
cat("MAF columns:\n"); print(names(maf))

#---------- ID handling ----------
sample_col <- NULL
if ("Tumor_Sample_Barcode" %in% names(maf)) {
  sample_col <- "Tumor_Sample_Barcode"
} else if ("submitter_id" %in% names(maf)) {
  sample_col <- "submitter_id"
} else {
  stop("No Tumor_Sample_Barcode or submitter_id column found in MAF.")
}

if (!"submitter_id" %in% names(maf)) {
  maf[, submitter_id := substr(get(sample_col), 1, 12)]
} else {
  # ensure 12-char
  maf[, submitter_id := substr(submitter_id, 1, 12)]
}

#---------- Basic columns ----------
required_cols <- c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "submitter_id")
missing <- setdiff(required_cols, names(maf))
if (length(missing) > 0) {
  stop("Missing required MAF columns: ", paste(missing, collapse = ", "))
}

# Keep SNVs only (1bp ref/alt)
maf_snvs <- maf[
  nchar(Reference_Allele) == 1 &
  nchar(Tumor_Seq_Allele2) == 1 &
  Reference_Allele %in% c("A","C","G","T") &
  Tumor_Seq_Allele2 %in% c("A","C","G","T")
]

cat("SNV rows retained:", nrow(maf_snvs), "\n")

#---------- Chromosome formatting ----------
maf_snvs[, Chromosome := as.character(Chromosome)]
maf_snvs[!grepl("^chr", Chromosome), Chromosome := paste0("chr", Chromosome)]

#---------- Get trinucleotide context ----------
hg <- BSgenome.Hsapiens.UCSC.hg38

get_tri <- function(chr, pos) {
  if (!chr %in% names(hg)) return(NA_character_)
  chr_len <- length(hg[[chr]])
  if (pos <= 1L || pos >= chr_len) return(NA_character_)
  as.character(subseq(hg[[chr]], start = pos - 1L, end = pos + 1L))
}

cat("Computing trinucleotide context...\n")
maf_snvs[, tri := get_tri(Chromosome, Start_Position), by = seq_len(nrow(maf_snvs))]
maf_snvs <- maf_snvs[!is.na(tri)]
cat("Rows with valid context:", nrow(maf_snvs), "\n")

#---------- Standardize to C or T ref (SBS96 convention) ----------
ref <- maf_snvs$Reference_Allele
alt <- maf_snvs$Tumor_Seq_Allele2
tri <- maf_snvs$tri

comp_vec <- c(A="T", C="G", G="C", T="A")

flip <- ref %in% c("A","G")

tri2 <- tri
if (any(flip)) {
  tri2[flip] <- as.character(reverseComplement(DNAStringSet(tri[flip])))
}

ref2 <- ref
ref2[flip] <- comp_vec[ref[flip]]

alt2 <- alt
alt2[flip] <- comp_vec[alt[flip]]

# sanity: only C/T now
keep_idx <- ref2 %in% c("C","T")
maf_snvs <- maf_snvs[keep_idx]
ref2 <- ref2[keep_idx]
alt2 <- alt2[keep_idx]
tri2 <- tri2[keep_idx]

# subtype = N[ref>alt]N
subtype <- paste0(substr(tri2, 1, 1), "[", ref2, ">", alt2, "]", substr(tri2, 3, 3))
maf_snvs[, subtype := subtype]
maf_snvs[, submitter_id := submitter_id]

cat("Unique SBS96 subtypes observed:", length(unique(subtype)), "\n")

#---------- Build SBS96 matrix ----------
tab <- maf_snvs[, .N, by = .(subtype, submitter_id)]
mat_dt <- dcast(tab, subtype ~ submitter_id, value.var = "N", fill = 0)
setorder(mat_dt, subtype)

cat("SBS96 matrix dimensions:", dim(mat_dt)[1], "x", dim(mat_dt)[2], "\n")
fwrite(mat_dt, out_mat, sep = "\t")
cat("Wrote SBS96 matrix to:\n  ", out_mat, "\n")

#---------- Signature-like scores per sample ----------
# total SNVs per sample
sample_cols <- setdiff(colnames(mat_dt), "subtype")
totals <- colSums(as.matrix(mat_dt[, ..sample_cols]))
total_dt <- data.table(submitter_id = sample_cols, total_snvs = as.integer(totals))

# Ageing-like: C>T at CpG (pattern N[C>T]G)
age_idx <- grepl("\\[C>T\\]", mat_dt$subtype) & grepl("G$", mat_dt$subtype)
age_mat <- as.matrix(mat_dt[age_idx, ..sample_cols])
age_counts <- colSums(age_mat)
age_prop <- ifelse(totals > 0, age_counts / totals, NA_real_)

# APOBEC-like: C>T or C>G at TpCpA/TpCpT (TpCpW)
apo_patterns <- c("T[C>T]A", "T[C>T]T", "T[C>G]A", "T[C>G]T")
apo_idx <- mat_dt$subtype %in% apo_patterns
apo_mat <- as.matrix(mat_dt[apo_idx, ..sample_cols])
apo_counts <- colSums(apo_mat)
apo_prop <- ifelse(totals > 0, apo_counts / totals, NA_real_)

sig_dt <- data.table(
  submitter_id        = sample_cols,
  total_snvs          = as.integer(totals),
  Ageing_CpG_prop     = as.numeric(age_prop),
  APOBEC_prop         = as.numeric(apo_prop)
)

cat("Signature score summary (first 5 rows):\n")
print(head(sig_dt))

fwrite(sig_dt, out_sc, sep = "\t")
cat("Wrote signature scores to:\n  ", out_sc, "\n")
cat("=== DONE SBS96 + Ageing/APOBEC score computation ===\n")
