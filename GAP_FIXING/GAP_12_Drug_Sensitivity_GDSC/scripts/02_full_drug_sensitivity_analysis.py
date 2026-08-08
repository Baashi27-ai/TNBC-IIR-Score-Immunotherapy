#!/usr/bin/env python3
"""
Script: 02_full_drug_sensitivity_analysis.py
Purpose: Correlate IIR score with ALL drug sensitivity in DepMap/PRISM
Author: Bhaskararao Ch (Baashi)
GAP: 12 — Drug Sensitivity Correlation (Full Analysis)
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')
from tqdm import tqdm

# =============================================================================
# 1. SETUP
# =============================================================================

os.chdir("D:/Baashi/TNBC_project/Immune_Spatial_Immunotherapy/GAP_FIXING/GAP_12_Drug_Sensitivity_GDSC")

print("=" * 60)
print(" FULL DRUG SENSITIVITY CORRELATION")
print(" IIR vs ALL Drug AUC — DepMap/PRISM")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

os.makedirs("results", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# =============================================================================
# 2. LOAD DATA
# =============================================================================

print("Loading data...")

# Load IIR scores
iir_file = "../GAP_02_ICB_Treated_Cohort_IMpassion/results/depmap/depmap_iir_scores.csv"

if not os.path.exists(iir_file):
    iir_file = "../../GAP_02_ICB_Treated_Cohort_IMpassion/results/depmap/depmap_iir_scores.csv"

if not os.path.exists(iir_file):
    print("ERROR: IIR file not found")
    exit(1)

iir_df = pd.read_csv(iir_file)
print(f"  IIR scores: {len(iir_df)} cell lines")

# Load PRISM drug sensitivity data
prism_file = "D:/Baashi/TNBC_project/Biomarker_Verification/inputs/external_cohorts/depmap/depmap_prism_auc.csv"

if not os.path.exists(prism_file):
    print("ERROR: PRISM file not found")
    exit(1)

prism_df = pd.read_csv(prism_file, low_memory=False)
print(f"  PRISM data: {prism_df.shape[0]} drugs × {prism_df.shape[1]-1} cell lines")

# Load cell line metadata
meta_file = "D:/Baashi/TNBC_project/Biomarker_Verification/inputs/external_cohorts/depmap/depmap_model.csv"

if os.path.exists(meta_file):
    meta_df = pd.read_csv(meta_file)
    print(f"  Metadata: {len(meta_df)} cell lines")
else:
    meta_df = None

# =============================================================================
# 3. MERGE IIR WITH DRUG DATA
# =============================================================================

print("\n" + "=" * 60)
print(" MERGING IIR WITH DRUG DATA")
print("=" * 60)

drug_col = prism_df.columns[0]
cell_cols = [col for col in prism_df.columns if col != drug_col]

iir_ids = iir_df['cell_line_id'].tolist()
prism_ids = cell_cols
overlap = set(iir_ids) & set(prism_ids)
overlap_list = list(overlap)
print(f"  Overlap: {len(overlap_list)} cell lines")

iir_filtered = iir_df[iir_df['cell_line_id'].isin(overlap_list)]
prism_filtered = prism_df[[drug_col] + overlap_list]

print(f"  IIR filtered: {len(iir_filtered)} cell lines")
print(f"  PRISM filtered: {prism_filtered.shape[0]} drugs × {prism_filtered.shape[1]-1} cell lines")

# =============================================================================
# 4. DRUG SENSITIVITY CORRELATIONS — ALL DRUGS
# =============================================================================

print("\n" + "=" * 60)
print(" DRUG SENSITIVITY CORRELATIONS — ALL DRUGS")
print(" IIR vs Drug AUC (lower AUC = higher sensitivity)")
print("=" * 60)

correlation_results = []
failed_count = 0

# Iterate through ALL drugs
total_drugs = len(prism_filtered)
print(f"\n  Processing {total_drugs} drugs... (this may take 30-60 minutes)")

for idx, row in tqdm(prism_filtered.iterrows(), total=total_drugs, desc="Processing drugs"):
    drug_name = row[drug_col]
    
    # Get AUC values - convert to numeric
    drug_auc = pd.to_numeric(row[overlap_list], errors='coerce').values
    
    # Remove NA
    valid_mask = ~np.isnan(drug_auc)
    valid_auc = drug_auc[valid_mask]
    valid_ids = np.array(overlap_list)[valid_mask]
    
    if len(valid_auc) < 10:
        failed_count += 1
        continue
    
    # Get IIR scores for valid cell lines
    iir_scores = []
    for cell_id in valid_ids:
        iir_val = iir_filtered[iir_filtered['cell_line_id'] == cell_id]['IIR_score_norm'].values
        if len(iir_val) > 0:
            iir_scores.append(iir_val[0])
        else:
            iir_scores.append(np.nan)
    
    iir_scores = np.array(iir_scores)
    valid_mask2 = ~np.isnan(iir_scores)
    
    if np.sum(valid_mask2) < 10:
        failed_count += 1
        continue
    
    valid_auc2 = valid_auc[valid_mask2]
    valid_iir = iir_scores[valid_mask2]
    
    # Spearman correlation
    try:
        corr, p = spearmanr(valid_iir, valid_auc2)
        correlation_results.append({
            'drug': drug_name,
            'correlation': corr,
            'p_value': p,
            'n': len(valid_iir)
        })
    except Exception as e:
        failed_count += 1
        continue

# Convert to DataFrame
corr_df = pd.DataFrame(correlation_results)
corr_df = corr_df.sort_values('correlation', ascending=True)

print(f"\n  Total correlations computed: {len(corr_df)}")
print(f"  Failed/Skipped: {failed_count}")

sig_corr = corr_df[corr_df['p_value'] < 0.05]
print(f"  Significant correlations (p < 0.05): {len(sig_corr)} ({len(sig_corr)/len(corr_df)*100:.1f}%)")

# FDR correction
from statsmodels.stats.multitest import multipletests
if len(corr_df) > 0:
    p_values = corr_df['p_value'].values
    reject, p_adjusted, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    corr_df['p_adjusted'] = p_adjusted
    corr_df['significant_fdr'] = reject
    sig_fdr = corr_df[corr_df['significant_fdr'] == True]
    print(f"  Significant after FDR correction: {len(sig_fdr)}")

# =============================================================================
# 5. TOP DRUGS
# =============================================================================

print("\n" + "=" * 60)
print(" TOP DRUGS")
print("=" * 60)

if len(corr_df) > 0:
    print("\n  Top 20 drugs with NEGATIVE correlation (IIR-high = more sensitive):")
    for _, row in corr_df.head(20).iterrows():
        fdr_status = "✅" if row.get('significant_fdr', False) else ""
        print(f"    {fdr_status} {row['drug']}: r={row['correlation']:.4f}, p={row['p_value']:.4f}, p_adj={row.get('p_adjusted', np.nan):.4f}, n={row['n']}")

# =============================================================================
# 6. IDENTIFY IMMUNE-MODULATING DRUGS
# =============================================================================

print("\n" + "=" * 60)
print(" IMMUNE-MODULATING DRUGS")
print("=" * 60)

immune_keywords = [
    'BRD', 'JQ1', 'IBET', 'BET', 'OTX', 'CPI', 'GSK',
    'PIK3', 'BYL', 'alpelisib', 'taselisib', 'buparlisib', 'BKM',
    'HDAC', 'vorinostat', 'panobinostat', 'romidepsin',
    'lenalidomide', 'pomalidomide',
    'TLR', 'STING', 'IFN', 'IL-2',
    'CXCR', 'CCR', 'TGFB', 'PD-1', 'PD-L1', 'CTLA-4'
]

immune_drugs = []

for _, row in corr_df.iterrows():
    drug_str = str(row['drug']).upper()
    for keyword in immune_keywords:
        if keyword.upper() in drug_str:
            immune_drugs.append(row)
            break

print(f"\n  Immune-modulating drugs in correlations: {len(immune_drugs)}")

if len(immune_drugs) > 0:
    immune_df = pd.DataFrame(immune_drugs)
    immune_df = immune_df.sort_values('correlation', ascending=True)
    print("\n  Immune-modulating drugs with negative correlation (IIR-high = more sensitive):")
    for _, row in immune_df.head(20).iterrows():
        print(f"    {row['drug']}: r={row['correlation']:.4f}, p={row['p_value']:.4f}")

# =============================================================================
# 7. DISTRIBUTION OF CORRELATIONS
# =============================================================================

print("\n" + "=" * 60)
print(" DISTRIBUTION OF CORRELATIONS")
print("=" * 60)

if len(corr_df) > 0:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Histogram
    ax = axes[0]
    ax.hist(corr_df['correlation'], bins=50, edgecolor='black', alpha=0.7)
    ax.axvline(0, color='red', linestyle='--')
    ax.set_xlabel('Correlation (Spearman)')
    ax.set_ylabel('Number of Drugs')
    ax.set_title('Distribution of IIR-Drug Correlations')
    ax.grid(True, alpha=0.3)
    
    # P-value distribution
    ax = axes[1]
    ax.hist(corr_df['p_value'], bins=50, edgecolor='black', alpha=0.7)
    ax.axvline(0.05, color='red', linestyle='--', label='p = 0.05')
    ax.set_xlabel('p-value')
    ax.set_ylabel('Number of Drugs')
    ax.set_title('Distribution of p-values')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    fig_path = "results/figures/correlation_distribution.png"
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"  ✓ Distribution plot saved: {fig_path}")

# =============================================================================
# 8. SAVE RESULTS
# =============================================================================

print("\nSaving results...")

corr_df.to_csv("results/drug_correlations_full.tsv", sep='\t', index=False)
print("  ✓ Full drug correlations saved")

if len(immune_drugs) > 0:
    immune_df.to_csv("results/immune_drug_correlations_full.tsv", sep='\t', index=False)
    print("  ✓ Immune drug correlations saved")

# Summary
summary = f"""FULL DRUG SENSITIVITY CORRELATION SUMMARY
========================================
Cell lines: {len(iir_filtered)}
Drugs analyzed: {len(corr_df)}
Significant correlations (p < 0.05): {len(sig_corr)} ({len(sig_corr)/len(corr_df)*100:.1f}%)
Significant after FDR: {len(sig_fdr) if 'sig_fdr' in locals() else 0}

Top 5 drugs (negative correlation):
"""
for _, row in corr_df.head(5).iterrows():
    summary += f"  {row['drug']}: r={row['correlation']:.4f}, p={row['p_value']:.4f}\n"

summary += f"\nImmune-modulating drugs: {len(immune_drugs)}\n"

with open("results/drug_sensitivity_full_summary.txt", "w") as f:
    f.write(summary)

print("  ✓ Summary saved")

# =============================================================================
# 9. SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print(" SUMMARY — FULL DRUG SENSITIVITY CORRELATION")
print("=" * 60)

print(f"\n  Cell lines: {len(iir_filtered)}")
print(f"  Drugs analyzed: {len(corr_df)}")
print(f"  Significant correlations (p < 0.05): {len(sig_corr)} ({len(sig_corr)/len(corr_df)*100:.1f}%)")
if 'sig_fdr' in locals():
    print(f"  Significant after FDR: {len(sig_fdr)}")

if len(corr_df) > 0:
    print(f"\n  Best correlation: {corr_df.iloc[0]['drug']}")
    print(f"    r = {corr_df.iloc[0]['correlation']:.4f}, p = {corr_df.iloc[0]['p_value']:.4f}")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ FULL DRUG SENSITIVITY CORRELATION COMPLETE")