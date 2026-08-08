#!/usr/bin/env python3
"""
Script: 05_compute_dpi.py
Purpose: Compute Drug Penetration Index (DPI) from ImmuneScore and spatial metrics
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
print(" DRUG PENETRATION INDEX (DPI) - METABRIC TNBC")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

# Create output directory
os.makedirs("results/iir", exist_ok=True)

# =============================================================================
# 2. LOAD INPUTS
# =============================================================================

print("Loading inputs...")

# Load ImmuneScore
immune_file = "results/immune_scores/immune_scores_metabric.tsv"
if not os.path.exists(immune_file):
    print(f"ERROR: ImmuneScore file not found: {immune_file}")
    print("Please run 02_compute_immune_scores.py first.")
    exit(1)

immune_df = pd.read_csv(immune_file, sep='\t')
print(f"  ImmuneScore: {immune_df.shape[0]} samples")

# Load Spatial Metrics
spatial_file = "results/spatial_metrics/spatial_metrics_metabric.tsv"
if not os.path.exists(spatial_file):
    print(f"ERROR: Spatial metrics file not found: {spatial_file}")
    print("Please run 04_compute_spatial_metrics.py first.")
    exit(1)

spatial_df = pd.read_csv(spatial_file, sep='\t')
print(f"  Spatial metrics: {spatial_df.shape[0]} samples")

# =============================================================================
# 3. MERGE DATA
# =============================================================================

# Merge on sample_id
merged = pd.merge(immune_df, spatial_df, on='sample_id')
print(f"  Merged: {merged.shape[0]} samples")

# =============================================================================
# 4. COMPUTE DPI
# =============================================================================

print("\nComputing Drug Penetration Index (DPI)...")

# DPI Formula (from Phase V):
# DPI = (ImmuneScore_norm × Immune_Stroma_ratio) / (1 + Immune_Exclusion_index)
# This captures both immune infiltration AND spatial accessibility

merged['DPI'] = (
    merged['ImmuneScore_norm'] * merged['Immune_Stroma_ratio']
) / (1 + merged['Immune_Exclusion_index'])

# Normalize DPI to 0-1 range
min_dpi = merged['DPI'].min()
max_dpi = merged['DPI'].max()
merged['DPI_norm'] = (merged['DPI'] - min_dpi) / (max_dpi - min_dpi + 1e-10)

print(f"  DPI range: {merged['DPI'].min():.4f} - {merged['DPI'].max():.4f}")
print(f"  DPI_norm range: {merged['DPI_norm'].min():.4f} - {merged['DPI_norm'].max():.4f}")

# =============================================================================
# 5. DPI TERTILES
# =============================================================================

print("\nClassifying DPI tertiles...")

tertiles = np.percentile(merged['DPI_norm'], [33.33, 66.67])
print(f"  Tertile thresholds: {tertiles[0]:.4f}, {tertiles[1]:.4f}")

merged['DPI_tertile'] = pd.cut(
    merged['DPI_norm'],
    bins=[-np.inf, tertiles[0], tertiles[1], np.inf],
    labels=['Low', 'Mid', 'High']
)

print(f"\n  DPI tertile distribution:")
print(merged['DPI_tertile'].value_counts().sort_index())

# =============================================================================
# 6. CREATE OUTPUT DATAFRAME
# =============================================================================

result_df = merged[[
    'sample_id',
    'ImmuneScore_norm',
    'Immune_Stroma_ratio',
    'Immune_Exclusion_index',
    'DPI',
    'DPI_norm',
    'DPI_tertile'
]].copy()

# =============================================================================
# 7. SAVE RESULTS
# =============================================================================

output_file = "results/iir/dpi_metabric.tsv"
result_df.to_csv(output_file, sep='\t', index=False)

print(f"\n  ✓ Saved to: {output_file}")

# =============================================================================
# 8. SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print(" SUMMARY")
print("=" * 60)

print(f"\nDPI summary:")
print(f"  Min: {result_df['DPI'].min():.4f}")
print(f"  Q1: {result_df['DPI'].quantile(0.25):.4f}")
print(f"  Median: {result_df['DPI'].median():.4f}")
print(f"  Q3: {result_df['DPI'].quantile(0.75):.4f}")
print(f"  Max: {result_df['DPI'].max():.4f}")

print(f"\nDPI_norm summary:")
print(f"  Min: {result_df['DPI_norm'].min():.4f}")
print(f"  Q1: {result_df['DPI_norm'].quantile(0.25):.4f}")
print(f"  Median: {result_df['DPI_norm'].median():.4f}")
print(f"  Q3: {result_df['DPI_norm'].quantile(0.75):.4f}")
print(f"  Max: {result_df['DPI_norm'].max():.4f}")

# Check correlations
print(f"\nCorrelations with DPI_norm:")
print(f"  vs ImmuneScore_norm: {result_df['DPI_norm'].corr(result_df['ImmuneScore_norm']):.4f}")
print(f"  vs Immune_Stroma_ratio: {result_df['DPI_norm'].corr(result_df['Immune_Stroma_ratio']):.4f}")
print(f"  vs Immune_Exclusion_index: {result_df['DPI_norm'].corr(result_df['Immune_Exclusion_index']):.4f}")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ DPI COMPUTATION COMPLETE")
print("   Next: 06_compute_iir.py")