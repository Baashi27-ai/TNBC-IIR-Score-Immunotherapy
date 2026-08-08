#!/usr/bin/env python3
"""
Script: 02_compute_immune_scores.py
Purpose: Compute ImmuneScore from xCell results
Author: Bhaskararao Ch (Baashi)
"""

import os
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# 1. SETUP
# =============================================================================

os.chdir("D:/Baashi/TNBC_project/Immune_Spatial_Immunotherapy/GAP_FIXING/GAP_01_External_Validation_METABRIC")

print("=" * 60)
print(" IMMUNESCORE COMPUTATION - METABRIC TNBC")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

# Create output directory
os.makedirs("results/immune_scores", exist_ok=True)

# =============================================================================
# 2. LOAD xCELL RESULTS
# =============================================================================

print("Loading xCell results...")

xcell_file = "results/xcell/xcell_scores_metabric.tsv"

if not os.path.exists(xcell_file):
    print(f"ERROR: xCell file not found: {xcell_file}")
    print("Please run 01_xcell_metabric.py first.")
    exit(1)

xcell_df = pd.read_csv(xcell_file, sep='\t')

print(f"  Loaded: {xcell_df.shape[0]} samples × {xcell_df.shape[1]} columns")

# Get sample IDs
sample_ids = xcell_df['sample_id'].tolist()

# =============================================================================
# 3. DEFINE IMMUNE CELL TYPES
# =============================================================================

# Core immune cell types for ImmuneScore (based on literature)
immune_cell_types = [
    'B cells',
    'CD4+ T-cells',
    'CD8+ T-cells',
    'T cells',
    'Tregs',
    'NK cells',
    'Macrophages M1',
    'Macrophages M2',
    'Monocytes',
    'Myeloid dendritic cells',
    'Plasmacytoid dendritic cells',
    'Eosinophils',
    'Mast cells',
    'Neutrophils'
]

# Check which are available
available_immune = [ct for ct in immune_cell_types if ct in xcell_df.columns]

print(f"\n  Immune cell types: {len(available_immune)}")
for ct in available_immune:
    print(f"    {ct}")

# =============================================================================
# 4. COMPUTE IMMUNESCORE
# =============================================================================

print("\nComputing ImmuneScore...")

# Method 1: Sum of all immune cell type scores
immune_scores = xcell_df[available_immune].sum(axis=1)

# Normalize to 0-1 scale
min_score = immune_scores.min()
max_score = immune_scores.max()
immune_scores_norm = (immune_scores - min_score) / (max_score - min_score + 1e-10)

# Method 2: Alternative - mean of z-scores (more robust)
# We'll use both and keep the sum-based one as main

print(f"  ImmuneScore range: {immune_scores.min():.4f} - {immune_scores.max():.4f}")
print(f"  Normalized range: {immune_scores_norm.min():.4f} - {immune_scores_norm.max():.4f}")

# =============================================================================
# 5. CLASSIFY HOT/INTERMEDIATE/COLD
# =============================================================================

print("\nClassifying tumors...")

# Use tertiles matching your TCGA approach
tertiles = np.percentile(immune_scores_norm, [33.33, 66.67])

print(f"  Tertile thresholds: {tertiles[0]:.4f}, {tertiles[1]:.4f}")

# Assign groups
immune_groups = []
for score in immune_scores_norm:
    if score <= tertiles[0]:
        immune_groups.append('Cold')
    elif score <= tertiles[1]:
        immune_groups.append('Intermediate')
    else:
        immune_groups.append('Hot')

# =============================================================================
# 6. CREATE OUTPUT DATAFRAME
# =============================================================================

result_df = pd.DataFrame({
    'sample_id': sample_ids,
    'ImmuneScore_raw': immune_scores,
    'ImmuneScore_norm': immune_scores_norm,
    'immune_group': immune_groups
})

# Add individual cell type scores
for ct in available_immune:
    result_df[ct] = xcell_df[ct].values

# =============================================================================
# 7. SAVE RESULTS
# =============================================================================

output_file = "results/immune_scores/immune_scores_metabric.tsv"
result_df.to_csv(output_file, sep='\t', index=False)

print(f"\n  ✓ Saved to: {output_file}")

# =============================================================================
# 8. SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print(" SUMMARY")
print("=" * 60)

print(f"\nImmuneScore summary:")
print(f"  Min: {result_df['ImmuneScore_norm'].min():.4f}")
print(f"  Q1: {result_df['ImmuneScore_norm'].quantile(0.25):.4f}")
print(f"  Median: {result_df['ImmuneScore_norm'].median():.4f}")
print(f"  Q3: {result_df['ImmuneScore_norm'].quantile(0.75):.4f}")
print(f"  Max: {result_df['ImmuneScore_norm'].max():.4f}")

print(f"\nImmune Group distribution:")
group_counts = result_df['immune_group'].value_counts()
for group in ['Hot', 'Intermediate', 'Cold']:
    count = group_counts.get(group, 0)
    pct = count / len(result_df) * 100
    print(f"  {group}: {count} ({pct:.1f}%)")

# Check correlation between cell types
print(f"\nTop correlations with ImmuneScore:")
immune_cols = [ct for ct in available_immune if ct in result_df.columns]
correlations = result_df[immune_cols].corrwith(result_df['ImmuneScore_norm'])
top_corrs = correlations.abs().sort_values(ascending=False).head(5)
for ct, corr in top_corrs.items():
    print(f"  {ct}: {corr:.3f}")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ IMMUNESCORE COMPUTATION COMPLETE")
print("   Next: 03_compute_pd1_signature.py")