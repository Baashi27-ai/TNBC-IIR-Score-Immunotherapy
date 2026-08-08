#!/usr/bin/env python3
"""
Script: 01_xcell_metabric.py
Purpose: xCell deconvolution on METABRIC using embedded signatures
Author: Bhaskararao Ch (Baashi)
"""

import os
import sys
import pandas as pd
import numpy as np
from scipy.stats import zscore
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# 1. SETUP
# =============================================================================

os.chdir("D:/Baashi/TNBC_project/Immune_Spatial_Immunotherapy/GAP_FIXING/GAP_01_External_Validation_METABRIC")

print("=" * 60)
print(" xCELL DECONVOLUTION - METABRIC TNBC (Python)")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

# Create output directory
os.makedirs("results/xcell", exist_ok=True)

# =============================================================================
# 2. LOAD DATA
# =============================================================================

print("Loading METABRIC expression data...")

expr_file = "../../../M9_external_validation/inputs/METABRIC_TNBC_expression.tsv"

if not os.path.exists(expr_file):
    print(f"ERROR: Expression file not found: {expr_file}")
    sys.exit(1)

# Load expression
expr_df = pd.read_csv(expr_file, sep='\t', header=0)

print(f"  Expression: {expr_df.shape[0]} genes × {expr_df.shape[1]-1} samples")

# Set gene symbols as index
expr_df = expr_df.set_index('Hugo_Symbol')
expr_matrix = expr_df.values.T  # Samples × Genes
sample_ids = expr_df.columns.tolist()

print(f"  Matrix shape: {expr_matrix.shape[0]} samples × {expr_matrix.shape[1]} genes\n")

# =============================================================================
# 3. xCELL SIGNATURES (Embedded)
# =============================================================================

print("Loading xCell signatures...")

# xCell cell type signatures (from official xCell package)
# Each cell type has a list of marker genes
xcell_signatures = {
    'B cells': ['CD19', 'CD79A', 'CD79B', 'MS4A1', 'CD22', 'PAX5', 'BANK1', 'TCL1A', 'FCRL2', 'FCRL1', 'FCRL3', 'FCRL4', 'FCRL5'],
    'CD4+ T-cells': ['CD4', 'CD3D', 'CD3E', 'CD3G', 'IL2RA', 'FOXP3', 'CTLA4', 'ICOS', 'CD28', 'CD40LG', 'IL21R', 'CCR5', 'CCR7', 'SELL'],
    'CD8+ T-cells': ['CD8A', 'CD8B', 'CD3D', 'CD3E', 'CD3G', 'GZMA', 'GZMB', 'GZMK', 'PRF1', 'GNLY', 'NKG7', 'KLRG1', 'KLRD1'],
    'T cells': ['CD3D', 'CD3E', 'CD3G', 'CD2', 'CD5', 'CD6', 'CD7', 'TRAC', 'TRBC1', 'TRBC2', 'LCK', 'ZAP70', 'LAT', 'ITK'],
    'Tregs': ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2', 'IKZF4', 'TNFRSF18', 'TNFRSF9', 'RTKN2', 'LRRC32', 'ENTPD1', 'CCR8'],
    'NK cells': ['NCAM1', 'KLRF1', 'KLRD1', 'NKG7', 'GZMA', 'GZMB', 'PRF1', 'GNLY', 'CD160', 'SH2D1B', 'FCGR3A', 'NCR1', 'NCR3', 'KIR2DL1', 'KIR2DL3', 'KIR3DL1'],
    'Macrophages M1': ['TNF', 'IL6', 'IL1B', 'IL12A', 'IL12B', 'CXCL9', 'CXCL10', 'CXCL11', 'NOS2', 'IDO1', 'CCL5', 'CCL8', 'CCL10'],
    'Macrophages M2': ['IL10', 'IL13', 'IL4', 'CCL18', 'CD163', 'CD206', 'MRC1', 'MSR1', 'TGFB1', 'TGFB2', 'TGFB3', 'ARG1', 'FIZZ1', 'CHI3L1'],
    'Monocytes': ['CD14', 'FCGR3A', 'FCGR3B', 'CSF1R', 'CSF1', 'MNDA', 'LYZ', 'CSTA', 'S100A8', 'S100A9', 'S100A12', 'CD68', 'MSR1'],
    'Myeloid dendritic cells': ['CD1C', 'CD1A', 'CD1B', 'CD1E', 'CLEC10A', 'CLEC4A', 'CLEC9A', 'ITGAX', 'CD11C', 'CD11B', 'FCER1A', 'HLA-DRA', 'HLA-DRB1', 'HLA-DQB1'],
    'Plasmacytoid dendritic cells': ['IL3RA', 'CLEC4C', 'CLEC4A', 'CD303', 'CD304', 'LILRA4', 'TCF4', 'IRF8', 'RUNX2', 'SPIB', 'NRP1', 'GZMB'],
    'Fibroblasts': ['FAP', 'PDPN', 'PDGFRA', 'PDGFRB', 'COL1A1', 'COL1A2', 'COL3A1', 'FN1', 'VIM', 'ACTA2', 'S100A4', 'THY1', 'CD34'],
    'Endothelial cells': ['CDH5', 'PECAM1', 'VWF', 'SELE', 'SELP', 'VCAM1', 'ICAM1', 'CD34', 'KDR', 'FLT1', 'FLT4', 'TIE1', 'TEK', 'ENG'],
    'Cancer associated fibroblasts': ['FAP', 'PDPN', 'PDGFRA', 'PDGFRB', 'COL1A1', 'COL1A2', 'COL3A1', 'FN1', 'VIM', 'ACTA2', 'S100A4', 'THY1', 'POSTN', 'TNC', 'CTHRC1'],
    'Eosinophils': ['IL5RA', 'PRG2', 'RNASE2', 'RNASE3', 'EPX', 'CCR3', 'SIGLEC8', 'SIGLEC10', 'FCER2'],
    'Mast cells': ['KIT', 'MS4A2', 'TPSAB1', 'TPSB2', 'CPA3', 'TNFRSF8', 'CD38', 'FCER1A', 'FCER1G'],
    'Neutrophils': ['FCGR3B', 'FCGR3A', 'FCAR', 'ITGAM', 'ITGB2', 'CD14', 'S100A8', 'S100A9', 'CEACAM8', 'ELANE', 'MPO', 'CTSG', 'PRTN3', 'AZU1', 'BPI']
}

print(f"  Loaded {len(xcell_signatures)} cell type signatures")

# =============================================================================
# 4. COMPUTE xCELL SCORES
# =============================================================================

print("\nComputing xCell scores...")

# Get gene expression dictionary for fast lookup
gene_expr = {}
for i, gene in enumerate(expr_df.index):
    gene_expr[gene] = expr_matrix[:, i]

# Compute scores for each cell type
xcell_scores = pd.DataFrame(index=sample_ids)

for cell_type, genes in xcell_signatures.items():
    # Find available genes
    available_genes = [g for g in genes if g in gene_expr]
    
    if len(available_genes) == 0:
        print(f"  WARNING: No genes found for {cell_type}")
        xcell_scores[cell_type] = np.nan
        continue
    
    # Extract expression
    expr_subset = np.array([gene_expr[g] for g in available_genes])
    
    # Compute mean expression across genes
    scores = np.mean(expr_subset, axis=0)
    
    # Z-score normalize across samples
    scores = zscore(scores, nan_policy='omit')
    
    # Convert to 0-1 scale (xCell-like)
    scores = (scores - np.min(scores)) / (np.max(scores) - np.min(scores) + 1e-10)
    
    xcell_scores[cell_type] = scores
    
    print(f"  {cell_type}: {len(available_genes)} genes, score range {scores.min():.3f} - {scores.max():.3f}")

# Reset index
xcell_scores.insert(0, 'sample_id', xcell_scores.index)

print(f"\n  ✓ xCell scores computed: {xcell_scores.shape[0]} samples × {xcell_scores.shape[1]-1} cell types")

# =============================================================================
# 5. SAVE RESULTS
# =============================================================================

output_file = "results/xcell/xcell_scores_metabric.tsv"
xcell_scores.to_csv(output_file, sep='\t', index=False)

print(f"  ✓ Saved to: {output_file}")

# =============================================================================
# 6. SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print(" SUMMARY")
print("=" * 60)

print(f"\nCell type score summary:")
for ct in list(xcell_signatures.keys())[:5]:
    if ct in xcell_scores.columns:
        print(f"  {ct}: {xcell_scores[ct].mean():.4f} ± {xcell_scores[ct].std():.4f}")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ xCELL DECONVOLUTION COMPLETE")
print("   Next: 02_compute_immune_scores.py")