#!/usr/bin/env python3
"""
Script: 03_compute_pd1_signature.py
Purpose: Compute PD1/PD-L1 signature from METABRIC expression
Author: Bhaskararao Ch (Baashi)
"""

import os
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
print(" PD1/PD-L1 SIGNATURE - METABRIC TNBC")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

# Create output directory
os.makedirs("results/immune_scores", exist_ok=True)

# =============================================================================
# 2. LOAD EXPRESSION DATA
# =============================================================================

print("Loading METABRIC expression data...")

expr_file = "../../../M9_external_validation/inputs/METABRIC_TNBC_expression.tsv"

if not os.path.exists(expr_file):
    print(f"ERROR: Expression file not found: {expr_file}")
    exit(1)

expr_df = pd.read_csv(expr_file, sep='\t', header=0)
expr_df = expr_df.set_index('Hugo_Symbol')

sample_ids = expr_df.columns.tolist()

print(f"  Loaded: {expr_df.shape[0]} genes × {expr_df.shape[1]} samples")

# =============================================================================
# 3. DEFINE PD1/PD-L1 SIGNATURE GENES
# =============================================================================

# PD1/PD-L1 pathway genes (matching your TCGA signature)
pd1_genes = [
    # Checkpoint genes
    'PDCD1', 'CD274', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT', 'BTLA',
    # Chemokines (T-cell recruitment)
    'CXCL9', 'CXCL10', 'CXCL11', 'CXCL13',
    # IFN signaling
    'IFNG', 'STAT1', 'IRF1',
    # CD8 T-cell markers
    'CD8A', 'CD8B', 'GZMA', 'GZMB', 'PRF1', 'GNLY', 'NKG7',
    # T-cell markers
    'CD3D', 'CD3E', 'CD3G', 'CD4',
    # Antigen presentation
    'HLA-DRA', 'HLA-DRB1', 'HLA-DQB1'
]

# Find which genes are available
available_genes = [g for g in pd1_genes if g in expr_df.index]

print(f"\n  PD1 signature genes: {len(pd1_genes)} total, {len(available_genes)} available")

if len(available_genes) < 10:
    print(f"  WARNING: Only {len(available_genes)} genes found. Results may be unstable.")

# =============================================================================
# 4. COMPUTE PD1 SIGNATURE
# =============================================================================

print("\nComputing PD1/PD-L1 signature...")

# Extract expression for available genes
pd1_expr = expr_df.loc[available_genes]

# Z-score normalize each gene across samples
# Use apply with axis=1 to normalize each row (gene)
pd1_z = pd1_expr.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

# Convert to DataFrame if it's a Series
if isinstance(pd1_z, pd.Series):
    pd1_z = pd1_z.to_frame().T

# Compute signature as mean of z-scores (matching your TCGA method)
pd1_signature = pd1_z.mean(axis=0)

# Create output dataframe
result_df = pd.DataFrame({
    'sample_id': sample_ids,
    'PD1_PDL1_signature': pd1_signature.values
})

# Add median split group (High/Low)
median_val = result_df['PD1_PDL1_signature'].median()
result_df['PD1_group'] = np.where(
    result_df['PD1_PDL1_signature'] > median_val,
    'High', 'Low'
)

print(f"  PD1 signature range: {result_df['PD1_PDL1_signature'].min():.4f} - {result_df['PD1_PDL1_signature'].max():.4f}")
print(f"  Median: {median_val:.4f}")
print(f"  High group: {(result_df['PD1_group'] == 'High').sum()} samples")
print(f"  Low group: {(result_df['PD1_group'] == 'Low').sum()} samples")

# =============================================================================
# 5. GENE-LEVEL CONTRIBUTIONS (for reference)
# =============================================================================

# Compute mean contribution of each gene
gene_contributions = pd1_z.mean(axis=1).sort_values(ascending=False)
print(f"\n  Top contributing genes:")
for gene, val in list(gene_contributions.head(5).items()):
    print(f"    {gene}: {val:.3f}")

# =============================================================================
# 6. SAVE RESULTS
# =============================================================================

output_file = "results/immune_scores/pd1_signature_metabric.tsv"
result_df.to_csv(output_file, sep='\t', index=False)

print(f"\n  ✓ Saved to: {output_file}")

# =============================================================================
# 7. SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print(" SUMMARY")
print("=" * 60)

print(f"\nPD1/PD-L1 signature summary:")
print(f"  Min: {result_df['PD1_PDL1_signature'].min():.4f}")
print(f"  Q1: {result_df['PD1_PDL1_signature'].quantile(0.25):.4f}")
print(f"  Median: {result_df['PD1_PDL1_signature'].median():.4f}")
print(f"  Q3: {result_df['PD1_PDL1_signature'].quantile(0.75):.4f}")
print(f"  Max: {result_df['PD1_PDL1_signature'].max():.4f}")

print(f"\nPD1 Group distribution:")
print(f"  High: {(result_df['PD1_group'] == 'High').sum()}")
print(f"  Low: {(result_df['PD1_group'] == 'Low').sum()}")

# Check correlation with ImmuneScore
print("\nChecking correlation with ImmuneScore...")
immune_file = "results/immune_scores/immune_scores_metabric.tsv"
if os.path.exists(immune_file):
    immune_df = pd.read_csv(immune_file, sep='\t')
    merged = pd.merge(result_df, immune_df, on='sample_id')
    corr = merged['PD1_PDL1_signature'].corr(merged['ImmuneScore_norm'])
    print(f"  PD1 signature vs ImmuneScore: r = {corr:.4f}")
else:
    print("  ImmuneScore file not found. Skipping correlation.")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ PD1/PD-L1 SIGNATURE COMPLETE")
print("   Next: 04_compute_spatial_metrics.py")