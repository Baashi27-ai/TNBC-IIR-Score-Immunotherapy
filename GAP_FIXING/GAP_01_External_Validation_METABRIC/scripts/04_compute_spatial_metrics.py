#!/usr/bin/env python3
"""
Script: 04_compute_spatial_metrics.py
Purpose: Compute spatial metrics from xCell results
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
print(" SPATIAL METRICS - METABRIC TNBC")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

# Create output directory
os.makedirs("results/spatial_metrics", exist_ok=True)

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
# 3. DEFINE CELL TYPES
# =============================================================================

# Immune cell types (from previous scripts)
immune_types = [
    'B cells', 'CD4+ T-cells', 'CD8+ T-cells', 'T cells', 'Tregs',
    'NK cells', 'Macrophages M1', 'Macrophages M2', 'Monocytes',
    'Myeloid dendritic cells', 'Plasmacytoid dendritic cells',
    'Eosinophils', 'Mast cells', 'Neutrophils'
]

# Stromal cell types
stromal_types = [
    'Fibroblasts', 'Endothelial cells', 'Cancer associated fibroblasts'
]

# Check which are available
available_immune = [ct for ct in immune_types if ct in xcell_df.columns]
available_stromal = [ct for ct in stromal_types if ct in xcell_df.columns]

print(f"\n  Immune types available: {len(available_immune)}")
print(f"  Stromal types available: {len(available_stromal)}")

# =============================================================================
# 4. COMPUTE SPATIAL METRICS
# =============================================================================

print("\nComputing spatial metrics...")

# 4.1 Immune Composite (sum of immune cell types)
if available_immune:
    immune_composite = xcell_df[available_immune].sum(axis=1)
else:
    immune_composite = pd.Series([0] * len(xcell_df))
    print("  WARNING: No immune cell types found")

# 4.2 Stroma Composite (sum of stromal cell types)
if available_stromal:
    stroma_composite = xcell_df[available_stromal].sum(axis=1)
else:
    stroma_composite = pd.Series([0] * len(xcell_df))
    print("  WARNING: No stromal cell types found")

print(f"  Immune composite range: {immune_composite.min():.4f} - {immune_composite.max():.4f}")
print(f"  Stroma composite range: {stroma_composite.min():.4f} - {stroma_composite.max():.4f}")

# 4.3 Immune-Stroma Ratio
total = immune_composite + stroma_composite
immune_stroma_ratio = np.where(
    total > 0,
    immune_composite / total,
    0.5  # Default when both are 0
)

print(f"  Immune-Stroma Ratio range: {immune_stroma_ratio.min():.4f} - {immune_stroma_ratio.max():.4f}")

# 4.4 Immune Exclusion Index (1 - Ratio)
immune_exclusion_index = 1 - immune_stroma_ratio
print(f"  Immune Exclusion Index range: {immune_exclusion_index.min():.4f} - {immune_exclusion_index.max():.4f}")

# =============================================================================
# 5. CREATE OUTPUT DATAFRAME
# =============================================================================

result_df = pd.DataFrame({
    'sample_id': sample_ids,
    'Immune_composite': immune_composite,
    'Stroma_composite': stroma_composite,
    'Immune_Stroma_ratio': immune_stroma_ratio,
    'Immune_Exclusion_index': immune_exclusion_index
})

# Add tertiles for each metric
for metric in ['Immune_Stroma_ratio', 'Immune_Exclusion_index']:
    tertiles = np.percentile(result_df[metric], [33.33, 66.67])
    labels = ['Low', 'Mid', 'High']
    result_df[f'{metric}_tertile'] = pd.cut(
        result_df[metric],
        bins=[-np.inf, tertiles[0], tertiles[1], np.inf],
        labels=labels
    )

# =============================================================================
# 6. SAVE RESULTS
# =============================================================================

output_file = "results/spatial_metrics/spatial_metrics_metabric.tsv"
result_df.to_csv(output_file, sep='\t', index=False)

print(f"\n  ✓ Saved to: {output_file}")

# =============================================================================
# 7. SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print(" SUMMARY")
print("=" * 60)

print(f"\nSpatial metrics summary:")
for metric in ['Immune_composite', 'Stroma_composite', 'Immune_Stroma_ratio', 'Immune_Exclusion_index']:
    print(f"  {metric}:")
    print(f"    Min: {result_df[metric].min():.4f}")
    print(f"    Median: {result_df[metric].median():.4f}")
    print(f"    Max: {result_df[metric].max():.4f}")

print(f"\nImmune_Stroma_ratio tertile distribution:")
print(result_df['Immune_Stroma_ratio_tertile'].value_counts().sort_index())

print(f"\nImmune_Exclusion_index tertile distribution:")
print(result_df['Immune_Exclusion_index_tertile'].value_counts().sort_index())

# Check correlation with ImmuneScore
print("\nChecking correlations with ImmuneScore...")
immune_file = "results/immune_scores/immune_scores_metabric.tsv"
if os.path.exists(immune_file):
    immune_df = pd.read_csv(immune_file, sep='\t')
    merged = pd.merge(result_df, immune_df, on='sample_id')
    
    corr_ratio = merged['Immune_Stroma_ratio'].corr(merged['ImmuneScore_norm'])
    corr_exclusion = merged['Immune_Exclusion_index'].corr(merged['ImmuneScore_norm'])
    
    print(f"  Immune_Stroma_ratio vs ImmuneScore: r = {corr_ratio:.4f}")
    print(f"  Immune_Exclusion_index vs ImmuneScore: r = {corr_exclusion:.4f}")
else:
    print("  ImmuneScore file not found. Skipping correlations.")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ SPATIAL METRICS COMPLETE")
print("   Next: 05_compute_dpi.py")