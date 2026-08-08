#!/usr/bin/env python3
"""
Script: 01_predict_neoantigens.py
Purpose: Predict neoantigens from MAF data (simplified approach)
Author: Bhaskararao Ch (Baashi)
GAP: 09 — Neoantigen Prediction
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# 1. SETUP
# =============================================================================

os.chdir("D:/Baashi/TNBC_project/Immune_Spatial_Immunotherapy/GAP_FIXING/GAP_09_Neoantigen_Prediction")

print("=" * 60)
print(" NEOANTIGEN PREDICTION")
print(" Simplified Approach — MAF-based")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

os.makedirs("results", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# =============================================================================
# 2. LOAD DATA
# =============================================================================

print("Loading data...")

# Load MAF
maf_file = "D:/Baashi/TNBC_project/data_raw/tcga/mutations/tcga_mutations_TNBCproxy.tsv"

if not os.path.exists(maf_file):
    print(f"ERROR: MAF file not found: {maf_file}")
    exit(1)

maf_df = pd.read_csv(maf_file, sep='\t', comment='#', low_memory=False)
print(f"  MAF: {len(maf_df)} variants")

# Load HLA escape data (contains HLA types, ImmuneScore, IIR)
escape_file = "../../Phase_X_HLA_ImmuneEscape/results/hla/HLA_escape_TNBC_like_TCGA.tsv"

if not os.path.exists(escape_file):
    escape_file = "../../../Immune_Spatial_Immunotherapy/Phase_X_HLA_ImmuneEscape/results/hla/HLA_escape_TNBC_like_TCGA.tsv"

if not os.path.exists(escape_file):
    print("ERROR: HLA escape file not found")
    exit(1)

escape_df = pd.read_csv(escape_file, sep='\t')
print(f"  HLA escape data: {len(escape_df)} samples")

# =============================================================================
# 3. EXTRACT HLA ALLELES
# =============================================================================

print("\n" + "=" * 60)
print(" EXTRACTING HLA ALLELES")
print("=" * 60)

# Check HLA columns in escape_df
hla_cols = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F']
available_hla = [c for c in hla_cols if c in escape_df.columns]
print(f"  HLA columns available: {available_hla}")

# Create HLA allele strings
def get_hla_alleles(row):
    alleles = []
    for col in available_hla:
        if pd.notna(row[col]) and row[col] != '':
            alleles.append(str(row[col]))
    return '|'.join(alleles) if alleles else 'Unknown'

escape_df['HLA_alleles'] = escape_df.apply(get_hla_alleles, axis=1)
print(f"  HLA alleles extracted")

# =============================================================================
# 4. SIMPLIFIED NEOANTIGEN PREDICTION
# =============================================================================

print("\n" + "=" * 60)
print(" SIMPLIFIED NEOANTIGEN PREDICTION")
print("=" * 60)

# Define variant types that can generate neoantigens
neo_variant_types = [
    'Missense_Mutation',
    'Frame_Shift_Del',
    'Frame_Shift_Ins',
    'Nonsense_Mutation',
    'Splice_Site',
    'In_Frame_Del',
    'In_Frame_Ins'
]

# Filter to neoantigen-relevant variants
neo_maf = maf_df[maf_df['Variant_Classification'].isin(neo_variant_types)]
print(f"  Neoantigen-relevant variants: {len(neo_maf)}")

# Count per sample
sample_col = 'Tumor_Sample_Barcode'
if sample_col not in neo_maf.columns:
    sample_col = [c for c in neo_maf.columns if 'sample' in c.lower()][0]

variant_counts = neo_maf[sample_col].value_counts()
print(f"  Samples with neoantigen-relevant variants: {len(variant_counts)}")

# Create neoantigen score (simplified: count of variants with neoantigen potential)
neo_df = pd.DataFrame({
    'sample_id': variant_counts.index,
    'neo_variant_count': variant_counts.values,
    'neo_score': np.log1p(variant_counts.values)  # log transform
})

print(f"  Neoantigen score range: {neo_df['neo_score'].min():.4f} - {neo_df['neo_score'].max():.4f}")

# =============================================================================
# 5. MERGE WITH IMMUNE DATA
# =============================================================================

print("\n" + "=" * 60)
print(" MERGING WITH IMMUNE DATA")
print("=" * 60)

# Extract patient ID from sample_id
def extract_patient_id(sample_id):
    parts = str(sample_id).split('-')
    if len(parts) >= 3:
        return '-'.join(parts[:3])
    return str(sample_id)

neo_df['patient_id'] = neo_df['sample_id'].apply(extract_patient_id)

# Merge with escape data
escape_df['patient_id'] = escape_df['submitter_id'].apply(extract_patient_id)

merged_df = pd.merge(neo_df, escape_df, on='patient_id', how='inner')
print(f"  Merged: {len(merged_df)} samples")

# =============================================================================
# 6. CORRELATIONS: NEOANTIGEN vs IMMUNESCORE and IIR
# =============================================================================

print("\n" + "=" * 60)
print(" CORRELATIONS: NEOANTIGEN vs IMMUNE")
print("=" * 60)

# Prepare data
corr_df = merged_df[['neo_score', 'neo_variant_count', 'ImmuneScore', 'IIR_score']].dropna()

# Correlations with neo_score
if len(corr_df) > 10:
    corr_immune, p_immune = spearmanr(corr_df['neo_score'], corr_df['ImmuneScore'])
    corr_iir, p_iir = spearmanr(corr_df['neo_score'], corr_df['IIR_score'])
    
    print(f"\n  Neoantigen score vs ImmuneScore:")
    print(f"    r = {corr_immune:.4f}, p = {p_immune:.4f}")
    
    print(f"\n  Neoantigen score vs IIR_score:")
    print(f"    r = {corr_iir:.4f}, p = {p_iir:.4f}")
    
    # Check success criteria
    if corr_immune > 0.4 and p_immune < 0.05:
        print(f"\n  ✅ Neoantigen vs ImmuneScore: PASS (r > 0.4, p < 0.05)")
    else:
        print(f"\n  ❌ Neoantigen vs ImmuneScore: FAIL (r = {corr_immune:.4f})")
    
    if corr_iir > 0.4 and p_iir < 0.05:
        print(f"  ✅ Neoantigen vs IIR: PASS (r > 0.4, p < 0.05)")
    else:
        print(f"  ❌ Neoantigen vs IIR: FAIL (r = {corr_iir:.4f})")

# =============================================================================
# 7. NEOANTIGEN vs APOBEC GROUPS
# =============================================================================

print("\n" + "=" * 60)
print(" NEOANTIGEN vs APOBEC GROUPS")
print("=" * 60)

if 'APOBEC_group' in merged_df.columns:
    apobec_high = merged_df[merged_df['APOBEC_group'] == 'APOBEC_high']
    apobec_low = merged_df[merged_df['APOBEC_group'] == 'APOBEC_low']
    
    if len(apobec_high) > 0 and len(apobec_low) > 0:
        median_high = apobec_high['neo_score'].median()
        median_low = apobec_low['neo_score'].median()
        
        print(f"\n  Neoantigen score by APOBEC group:")
        print(f"    APOBEC_high: median = {median_high:.4f}, n = {len(apobec_high)}")
        print(f"    APOBEC_low: median = {median_low:.4f}, n = {len(apobec_low)}")
        
        # Mann-Whitney U test
        stat, p = mannwhitneyu(apobec_high['neo_score'], apobec_low['neo_score'])
        print(f"  Mann-Whitney U test: p = {p:.4f}")

# =============================================================================
# 8. NEOANTIGEN vs IMMUNE SUBTYPES
# =============================================================================

print("\n" + "=" * 60)
print(" NEOANTIGEN vs IMMUNE SUBTYPES")
print("=" * 60)

if 'immune_subtype' in merged_df.columns:
    print("\n  Neoantigen score by immune subtype:")
    for subtype in ['BLIS', 'IM', 'BLIA']:
        subset = merged_df[merged_df['immune_subtype'] == subtype]
        if len(subset) > 0:
            median_val = subset['neo_score'].median()
            print(f"    {subtype}: median = {median_val:.4f}, n = {len(subset)}")

# =============================================================================
# 9. VISUALIZE
# =============================================================================

print("\n" + "=" * 60)
print(" VISUALIZING NEOANTIGEN CORRELATIONS")
print("=" * 60)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Plot 1: Neo vs ImmuneScore
ax = axes[0]
ax.scatter(corr_df['neo_score'], corr_df['ImmuneScore'], alpha=0.5)
ax.set_xlabel('Neoantigen Score')
ax.set_ylabel('ImmuneScore')
ax.set_title(f'Neoantigen vs ImmuneScore\nr = {corr_immune:.4f}, p = {p_immune:.4f}')
ax.grid(True, alpha=0.3)

# Plot 2: Neo vs IIR
ax = axes[1]
ax.scatter(corr_df['neo_score'], corr_df['IIR_score'], alpha=0.5)
ax.set_xlabel('Neoantigen Score')
ax.set_ylabel('IIR Score')
ax.set_title(f'Neoantigen vs IIR\nr = {corr_iir:.4f}, p = {p_iir:.4f}')
ax.grid(True, alpha=0.3)

plt.tight_layout()
fig_path = "results/figures/neoantigen_correlations.png"
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
print(f"  ✓ Correlation plots saved: {fig_path}")

# =============================================================================
# 10. SAVE RESULTS
# =============================================================================

print("\nSaving results...")

# Save merged data
merged_df.to_csv("results/neoantigen_data.tsv", sep='\t', index=False)
print(f"  ✓ Data saved: results/neoantigen_data.tsv")

# Create summary
summary = f"""NEOANTIGEN PREDICTION SUMMARY
========================================
Dataset: TCGA TNBC-like
Samples: {len(merged_df)}
Neoantigen-relevant variants: {len(neo_maf)}
Neoantigen score range: {neo_df['neo_score'].min():.4f} - {neo_df['neo_score'].max():.4f}

Correlations (Spearman):
  Neoantigen vs ImmuneScore: r = {corr_immune:.4f}, p = {p_immune:.4f}
  Neoantigen vs IIR: r = {corr_iir:.4f}, p = {p_iir:.4f}

Neoantigen by APOBEC group:
  APOBEC_high: median = {median_high:.4f}, n = {len(apobec_high)}
  APOBEC_low: median = {median_low:.4f}, n = {len(apobec_low)}
  Mann-Whitney p = {p:.4f}
"""

with open("results/neoantigen_summary.txt", "w") as f:
    f.write(summary)

print(f"  ✓ Summary saved: results/neoantigen_summary.txt")

# =============================================================================
# 11. SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print(" SUMMARY — NEOANTIGEN PREDICTION")
print("=" * 60)

print(f"\n  Samples: {len(merged_df)}")
print(f"  Neoantigen-relevant variants: {len(neo_maf)}")
print(f"  Neoantigen score range: {neo_df['neo_score'].min():.4f} - {neo_df['neo_score'].max():.4f}")

print(f"\n  Correlations:")
print(f"    Neoantigen vs ImmuneScore: r = {corr_immune:.4f}, p = {p_immune:.4f}")
print(f"    Neoantigen vs IIR: r = {corr_iir:.4f}, p = {p_iir:.4f}")

if corr_immune > 0.4 and p_immune < 0.05:
    print("\n  ✅ Neoantigen vs ImmuneScore: PASS")
else:
    print(f"\n  ❌ Neoantigen vs ImmuneScore: FAIL (r = {corr_immune:.4f})")

if corr_iir > 0.4 and p_iir < 0.05:
    print("  ✅ Neoantigen vs IIR: PASS")
else:
    print(f"  ❌ Neoantigen vs IIR: FAIL (r = {corr_iir:.4f})")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ NEOANTIGEN PREDICTION COMPLETE")