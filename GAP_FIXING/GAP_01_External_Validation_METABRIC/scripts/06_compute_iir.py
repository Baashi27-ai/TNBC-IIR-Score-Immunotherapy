#!/usr/bin/env python3
"""
Script: 06_compute_iir.py
Purpose: Compute Integrated Immunotherapy Readiness (IIR) Score
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
print(" IIR SCORE - METABRIC TNBC")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

# Create output directory
os.makedirs("results/iir", exist_ok=True)

# =============================================================================
# 2. LOAD INPUTS
# =============================================================================

print("Loading inputs...")

# Load ImmuneScore
immune_df = pd.read_csv("results/immune_scores/immune_scores_metabric.tsv", sep='\t')
print(f"  ImmuneScore: {immune_df.shape[0]} samples")

# Load PD1 signature
pd1_df = pd.read_csv("results/immune_scores/pd1_signature_metabric.tsv", sep='\t')
print(f"  PD1 signature: {pd1_df.shape[0]} samples")

# Load Spatial metrics
spatial_df = pd.read_csv("results/spatial_metrics/spatial_metrics_metabric.tsv", sep='\t')
print(f"  Spatial metrics: {spatial_df.shape[0]} samples")

# Load DPI
dpi_df = pd.read_csv("results/iir/dpi_metabric.tsv", sep='\t')
print(f"  DPI: {dpi_df.shape[0]} samples")

# =============================================================================
# 3. MERGE DATA
# =============================================================================

print("\nMerging data...")

# Start with immune_df
merged = immune_df[['sample_id', 'ImmuneScore_norm']].copy()

# Add PD1 signature
merged = pd.merge(merged, pd1_df[['sample_id', 'PD1_PDL1_signature']], on='sample_id')

# Add spatial metrics
merged = pd.merge(merged, spatial_df[['sample_id', 'Immune_Exclusion_index']], on='sample_id')

# Add DPI
merged = pd.merge(merged, dpi_df[['sample_id', 'DPI_norm']], on='sample_id')

print(f"  Merged: {merged.shape[0]} samples")

# =============================================================================
# 4. COMPUTE IIR SCORE
# =============================================================================

print("\nComputing IIR Score...")

# 4.1 Z-score normalize each component
# ImmuneScore_norm
merged['ImmuneScore_z'] = zscore(merged['ImmuneScore_norm'])

# PD1_PDL1_signature
merged['PD1_z'] = zscore(merged['PD1_PDL1_signature'])

# -Immune_Exclusion_index (inverse because high exclusion = bad)
merged['Exclusion_inv'] = -merged['Immune_Exclusion_index']
merged['Exclusion_z'] = zscore(merged['Exclusion_inv'])

# DPI_norm
merged['DPI_z'] = zscore(merged['DPI_norm'])

# TMB and APOBEC are NOT available in METABRIC
# We'll set them to 0 (neutral) and note this limitation
merged['TMB_z'] = 0
merged['APOBEC_z'] = 0

print(f"  Components z-scored:")
print(f"    ImmuneScore_z: {merged['ImmuneScore_z'].min():.4f} - {merged['ImmuneScore_z'].max():.4f}")
print(f"    PD1_z: {merged['PD1_z'].min():.4f} - {merged['PD1_z'].max():.4f}")
print(f"    Exclusion_z: {merged['Exclusion_z'].min():.4f} - {merged['Exclusion_z'].max():.4f}")
print(f"    DPI_z: {merged['DPI_z'].min():.4f} - {merged['DPI_z'].max():.4f}")
print(f"    TMB_z: 0 (not available)")
print(f"    APOBEC_z: 0 (not available)")

# 4.2 Compute IIR as mean of z-scores
# Using 6 components (2 are placeholder zeros)
merged['IIR_score'] = (
    merged['ImmuneScore_z'] +
    merged['PD1_z'] +
    merged['Exclusion_z'] +
    merged['DPI_z'] +
    merged['TMB_z'] +
    merged['APOBEC_z']
) / 6

print(f"\n  IIR_score range: {merged['IIR_score'].min():.4f} - {merged['IIR_score'].max():.4f}")
print(f"  IIR_score mean: {merged['IIR_score'].mean():.4f}")
print(f"  IIR_score std: {merged['IIR_score'].std():.4f}")

# =============================================================================
# 5. IIR GROUPS (High/Intermediate/Poor)
# =============================================================================

print("\nClassifying IIR groups...")

# Use tertiles matching your TCGA approach
tertiles = np.percentile(merged['IIR_score'], [33.33, 66.67])
print(f"  Tertile thresholds: {tertiles[0]:.4f}, {tertiles[1]:.4f}")

merged['IIR_group'] = pd.cut(
    merged['IIR_score'],
    bins=[-np.inf, tertiles[0], tertiles[1], np.inf],
    labels=['Poor_ICB_ready', 'Intermediate', 'High_ICB_ready']
)

print(f"\n  IIR group distribution:")
print(merged['IIR_group'].value_counts().sort_index())

# =============================================================================
# 6. CREATE OUTPUT DATAFRAME
# =============================================================================

result_df = merged[[
    'sample_id',
    'ImmuneScore_norm',
    'PD1_PDL1_signature',
    'Immune_Exclusion_index',
    'DPI_norm',
    'ImmuneScore_z',
    'PD1_z',
    'Exclusion_z',
    'DPI_z',
    'TMB_z',
    'APOBEC_z',
    'IIR_score',
    'IIR_group'
]].copy()

# =============================================================================
# 7. SAVE RESULTS
# =============================================================================

output_file = "results/iir/iir_score_metabric.tsv"
result_df.to_csv(output_file, sep='\t', index=False)

print(f"\n  ✓ Saved to: {output_file}")

# =============================================================================
# 8. SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print(" SUMMARY")
print("=" * 60)

print(f"\nIIR_score summary:")
print(f"  Min: {result_df['IIR_score'].min():.4f}")
print(f"  Q1: {result_df['IIR_score'].quantile(0.25):.4f}")
print(f"  Median: {result_df['IIR_score'].median():.4f}")
print(f"  Q3: {result_df['IIR_score'].quantile(0.75):.4f}")
print(f"  Max: {result_df['IIR_score'].max():.4f}")

print(f"\nIIR group distribution:")
print(result_df['IIR_group'].value_counts().sort_index())

# Check correlations between components
print(f"\nCorrelations with IIR_score:")
components = ['ImmuneScore_z', 'PD1_z', 'Exclusion_z', 'DPI_z']
for comp in components:
    corr = result_df['IIR_score'].corr(result_df[comp])
    print(f"  {comp}: {corr:.4f}")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ IIR SCORE COMPLETE")
print("   Next: 07_validate_iir.py")