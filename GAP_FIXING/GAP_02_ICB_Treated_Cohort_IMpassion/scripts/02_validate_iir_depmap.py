#!/usr/bin/env python3
"""
Script: 02_validate_iir_depmap.py
Purpose: Validate IIR score against drug sensitivity in DepMap/PRISM
Author: Bhaskararao Ch (Baashi)
GAP: 02 — ICB-Treated Cohort Validation (Mechanistic)
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# 1. SETUP
# =============================================================================

os.chdir("D:/Baashi/TNBC_project/Immune_Spatial_Immunotherapy/GAP_FIXING/GAP_02_ICB_Treated_Cohort_IMpassion")

print("=" * 60)
print(" DEPMAP DRUG SENSITIVITY VALIDATION")
print(" IIR Score vs Drug Response")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

os.makedirs("results/depmap", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# =============================================================================
# 2. LOAD IIR SCORES
# =============================================================================

print("Loading IIR scores...")

iir_file = "results/depmap/depmap_iir_scores.csv"

if not os.path.exists(iir_file):
    print(f"ERROR: IIR file not found: {iir_file}")
    print("Please run 01_compute_iir_depmap.py first.")
    exit(1)

iir_df = pd.read_csv(iir_file)
print(f"  Loaded: {iir_df.shape[0]} cell lines")

# =============================================================================
# 3. LOAD PRISM DRUG SENSITIVITY DATA
# =============================================================================

print("\nLoading PRISM drug sensitivity data...")

prism_file = "D:/Baashi/TNBC_project/Biomarker_Verification/inputs/external_cohorts/depmap/depmap_prism_auc.csv"

if not os.path.exists(prism_file):
    print(f"ERROR: PRISM file not found: {prism_file}")
    exit(1)

prism_df = pd.read_csv(prism_file, low_memory=False)
print(f"  Loaded: {prism_df.shape[0]} rows × {prism_df.shape[1]} columns")

# =============================================================================
# 4. UNDERSTAND PRISM DATA STRUCTURE
# =============================================================================

print("\nUnderstanding PRISM data structure...")

# First column is drug name/identifier
drug_col = prism_df.columns[0]
print(f"  Drug column: {drug_col}")

# Other columns are cell line IDs (ACH-xxxxx)
cell_cols = [col for col in prism_df.columns if col != drug_col]
print(f"  Cell line columns: {len(cell_cols)}")

# Get drug names
drug_names = prism_df[drug_col].tolist()
print(f"  Number of drugs: {len(drug_names)}")
print(f"  First 5 drugs: {drug_names[:5]}")

# =============================================================================
# 5. FILTER TO DRUGS OF INTEREST
# =============================================================================

print("\nFiltering to drugs of interest...")

# Immunotherapy-related keywords
drug_keywords = [
    # Checkpoint inhibitors
    'PD-1', 'PD-L1', 'CTLA-4', 'PD1', 'PDL1', 'pembrolizumab', 'nivolumab', 
    'atezolizumab', 'avelumab', 'durvalumab', 'ipilimumab', 'tremelimumab',
    # BRD4 inhibitors (immune-modulating)
    'JQ1', 'IBET', 'BET', 'BRD4', 'OTX015', 'CPI-0610', 'GSK525762',
    # PI3K inhibitors
    'BYL719', 'alpelisib', 'taselisib', 'GDC-0032', 'buparlisib', 'BKM120',
    # Other immune modulators
    'lenalidomide', 'pomalidomide', 'thalidomide', 'TLR', 'STING', 'IFN',
    'IL-2', 'IL-15', 'CXCR', 'CCR', 'TGFB', 'TGF-β',
    # Epigenetic modifiers
    'HDAC', 'vorinostat', 'panobinostat', 'romidepsin', 'entinostat'
]

# Find drugs matching keywords
matching_drugs = []
for drug in drug_names:
    drug_str = str(drug)
    for keyword in drug_keywords:
        if keyword.lower() in drug_str.lower():
            matching_drugs.append(drug)
            break

print(f"  Found {len(matching_drugs)} drugs of interest")
for drug in matching_drugs[:10]:
    print(f"    {drug}")

if not matching_drugs:
    print("  No matching drugs found. Using all drugs with variance.")
    # Use all drugs
    matching_drugs = drug_names

# =============================================================================
# 6. EXTRACT DRUG DATA
# =============================================================================

print("\nExtracting drug sensitivity data...")

# Filter to matching drugs
drug_df = prism_df[prism_df[drug_col].isin(matching_drugs)]

if len(drug_df) == 0:
    print("  No data for matching drugs.")
    exit(1)

print(f"  Selected: {len(drug_df)} drugs")

# =============================================================================
# 7. TRANSPOSE TO CELL LINES × DRUGS
# =============================================================================

print("\nTransposing to cell lines × drugs format...")

# Set drug names as index
drug_df = drug_df.set_index(drug_col)

# Transpose: cell lines as rows, drugs as columns
drug_matrix = drug_df.T

print(f"  Transposed: {drug_matrix.shape[0]} cell lines × {drug_matrix.shape[1]} drugs")

# =============================================================================
# 8. MERGE IIR WITH DRUG DATA
# =============================================================================

print("\nMerging IIR scores with drug data...")

# Get cell line IDs from IIR
iir_ids = iir_df['cell_line_id'].tolist()
# Get cell line IDs from drug matrix
drug_ids = drug_matrix.index.tolist()

# Find overlap
overlap = set(iir_ids) & set(drug_ids)
print(f"  Overlap: {len(overlap)} cell lines")

if len(overlap) < 10:
    print("  WARNING: Limited overlap. Checking ID formats...")
    print(f"  IIR example IDs: {iir_ids[:5]}")
    print(f"  Drug matrix example IDs: {drug_ids[:5]}")

# Filter both datasets to overlap
iir_filtered = iir_df[iir_df['cell_line_id'].isin(overlap)]
drug_filtered = drug_matrix.loc[list(overlap)]

print(f"  IIR filtered: {len(iir_filtered)} cell lines")
print(f"  Drug filtered: {len(drug_filtered)} cell lines")

# Merge
merged_df = iir_filtered.copy()
for drug_col in drug_filtered.columns:
    merged_df[drug_col] = drug_filtered[drug_col].values

print(f"  Merged: {merged_df.shape[0]} cell lines × {merged_df.shape[1]} features")

# =============================================================================
# 9. DRUG SENSITIVITY CORRELATIONS
# =============================================================================

print("\n" + "=" * 60)
print(" DRUG SENSITIVITY CORRELATIONS")
print(" IIR_score vs Drug AUC (lower AUC = higher sensitivity)")
print("=" * 60)

correlation_results = []

# Get drug columns
drug_cols = [col for col in merged_df.columns if col not in ['cell_line_id', 'IIR_score_norm', 'IIR_group']]

print(f"\n  Analyzing {len(drug_cols)} drugs...")

for drug_col in drug_cols:
    try:
        # Filter out NA values
        data = merged_df[['IIR_score_norm', drug_col]].dropna()
        if len(data) < 10:
            continue
        
        # Spearman correlation
        corr, p_value = spearmanr(data['IIR_score_norm'], data[drug_col])
        
        # For AUC, negative correlation = higher IIR = more sensitive
        correlation_results.append({
            'drug': drug_col,
            'correlation': corr,
            'p_value': p_value,
            'n': len(data),
            'significance': 'p < 0.05' if p_value < 0.05 else 'NS'
        })
    except Exception as e:
        continue

# Sort by correlation
corr_df = pd.DataFrame(correlation_results)
if len(corr_df) > 0:
    corr_df = corr_df.sort_values('correlation', ascending=True)
    
    # Highlight significant correlations
    sig_corr = corr_df[corr_df['p_value'] < 0.05]
    print(f"\n  Significant correlations: {len(sig_corr)}/{len(corr_df)}")
    
    print(f"\n  Top 10 drugs with NEGATIVE correlation (IIR-high = more sensitive):")
    print(corr_df.head(10).to_string(index=False))
    
    print(f"\n  Top 10 drugs with POSITIVE correlation (IIR-high = less sensitive):")
    print(corr_df.tail(10).to_string(index=False))
    
    # Check if any immune-related drugs are in the top list
    immune_keywords = ['PD', 'PD1', 'PDL1', 'CTLA', 'JQ1', 'IBET', 'BRD4', 'BET']
    immune_drugs_in_top = []
    for _, row in corr_df.head(20).iterrows():
        for kw in immune_keywords:
            if kw.lower() in str(row['drug']).lower():
                immune_drugs_in_top.append(row['drug'])
                break
    
    if immune_drugs_in_top:
        print(f"\n  Immune-related drugs in top 20 (most sensitive):")
        for drug in immune_drugs_in_top[:10]:
            print(f"    {drug}")
    else:
        print("\n  No immune-related drugs found in top 20.")
        print("  This is expected — PRISM has limited immunotherapy drugs.")
        print("  BRD4/PI3K inhibitors are the best proxies for immune modulation.")
else:
    print("  No drug correlations computed.")

# =============================================================================
# 10. SAVE RESULTS
# =============================================================================

print("\nSaving results...")

if len(corr_df) > 0:
    corr_output = "results/depmap/depmap_drug_correlations.csv"
    corr_df.to_csv(corr_output, index=False)
    print(f"  ✓ Drug correlations saved: {corr_output}")

# Save merged data for further analysis
merged_output = "results/depmap/depmap_merged_data.csv"
merged_df.to_csv(merged_output, index=False)
print(f"  ✓ Merged data saved: {merged_output}")

# =============================================================================
# 11. SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print(" SUMMARY — DEPMAP DRUG VALIDATION")
print("=" * 60)

print(f"\n  Cell lines analyzed: {len(merged_df)}")
print(f"  Drugs analyzed: {len(drug_cols)}")

if len(corr_df) > 0:
    significant = corr_df[corr_df['p_value'] < 0.05]
    print(f"  Significant correlations: {len(significant)}")
    
    # Find the best negative correlation (more sensitive)
    best_neg = corr_df.loc[corr_df['correlation'].idxmin()]
    print(f"\n  Best negative correlation (IIR-high = more sensitive):")
    print(f"    Drug: {best_neg['drug']}")
    print(f"    Correlation: {best_neg['correlation']:.4f}")
    print(f"    p-value: {best_neg['p_value']:.4f}")
    print(f"    n: {best_neg['n']}")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ DEPMAP DRUG VALIDATION COMPLETE")