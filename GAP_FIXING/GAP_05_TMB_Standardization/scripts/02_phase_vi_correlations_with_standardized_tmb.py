#!/usr/bin/env python3
"""
Script: 02_phase_vi_correlations_with_standardized_tmb.py
Purpose: Re-run Phase VI correlations using standardized TMB
Author: Bhaskararao Ch (Baashi)
GAP: 05 — TMB Standardization
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# 1. SETUP
# =============================================================================

os.chdir("D:/Baashi/TNBC_project/Immune_Spatial_Immunotherapy/GAP_FIXING/GAP_05_TMB_Standardization")

print("=" * 60)
print(" PHASE VI CORRELATIONS — STANDARDIZED TMB")
print(" TMB vs Immune, Subtypes, Survival")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

os.makedirs("results/phase_vi", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# =============================================================================
# 2. LOAD DATA
# =============================================================================

print("Loading data...")

# Load standardized TMB
tmb_file = "results/tmb_standardized.tsv"
if not os.path.exists(tmb_file):
    print(f"ERROR: TMB file not found: {tmb_file}")
    exit(1)

tmb_df = pd.read_csv(tmb_file, sep='\t')
print(f"  Standardized TMB: {len(tmb_df)} samples")
print(f"  TMB median: {tmb_df['TMB_mut_per_Mb'].median():.2f} mut/Mb")

# Load HLA escape data (contains immune scores, subtypes, IIR)
escape_file = "../../Phase_X_HLA_ImmuneEscape/results/hla/HLA_escape_TNBC_like_TCGA.tsv"

if not os.path.exists(escape_file):
    escape_file = "../../../Immune_Spatial_Immunotherapy/Phase_X_HLA_ImmuneEscape/results/hla/HLA_escape_TNBC_like_TCGA.tsv"

if not os.path.exists(escape_file):
    print("ERROR: HLA escape file not found")
    exit(1)

escape_df = pd.read_csv(escape_file, sep='\t')
print(f"  HLA escape data: {len(escape_df)} samples")

# =============================================================================
# 3. MERGE DATA
# =============================================================================

print("\nMerging data...")

# Extract sample ID from TMB (Tumor_Sample_Barcode -> TCGA-XX-XXXX)
def extract_patient_id(sample_id):
    parts = str(sample_id).split('-')
    if len(parts) >= 3:
        return '-'.join(parts[:3])
    return str(sample_id)

tmb_df['patient_id'] = tmb_df['sample_id'].apply(extract_patient_id)

# Merge with escape data (using submitter_id as patient ID)
escape_df['patient_id'] = escape_df['submitter_id'].apply(extract_patient_id)

# Merge
merged_df = pd.merge(tmb_df, escape_df, on='patient_id', how='inner')
print(f"  Merged: {len(merged_df)} samples")

# =============================================================================
# 4. TMB vs ImmuneScore
# =============================================================================

print("\n" + "=" * 60)
print(" TMB vs IMMUNESCORE")
print("=" * 60)

if 'ImmuneScore' in merged_df.columns:
    corr, p = spearmanr(merged_df['TMB_mut_per_Mb'], merged_df['ImmuneScore'])
    print(f"\n  Spearman correlation: r = {corr:.4f}, p = {p:.4f}")
    
    # Scatter plot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(merged_df['TMB_mut_per_Mb'], merged_df['ImmuneScore'], alpha=0.5)
    ax.set_xlabel('TMB (mut/Mb)')
    ax.set_ylabel('ImmuneScore')
    ax.set_title(f'TMB vs ImmuneScore\nr = {corr:.4f}, p = {p:.4f}')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    fig_path = "results/figures/tmb_vs_immunescore.png"
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"  ✓ Plot saved: {fig_path}")

# =============================================================================
# 5. TMB vs Immune Subtypes
# =============================================================================

print("\n" + "=" * 60)
print(" TMB vs IMMUNE SUBTYPES")
print("=" * 60)

if 'immune_subtype' in merged_df.columns:
    # TMB by subtype
    print("\n  TMB by immune subtype:")
    subtype_data = {}
    for subtype in ['BLIS', 'IM', 'BLIA']:
        subset = merged_df[merged_df['immune_subtype'] == subtype]
        if len(subset) > 0:
            mean_val = subset['TMB_mut_per_Mb'].mean()
            std_val = subset['TMB_mut_per_Mb'].std()
            median_val = subset['TMB_mut_per_Mb'].median()
            print(f"    {subtype}: mean={mean_val:.2f} ± {std_val:.2f}, median={median_val:.2f}, n={len(subset)}")
            subtype_data[subtype] = subset['TMB_mut_per_Mb'].values
    
    # Boxplot (fixed: use tick labels instead of labels parameter)
    fig, ax = plt.subplots(figsize=(8, 6))
    data = [subtype_data.get(st, []) for st in ['BLIS', 'IM', 'BLIA']]
    bp = ax.boxplot(data, patch_artist=True)
    ax.set_xticklabels(['BLIS', 'IM', 'BLIA'])
    ax.set_ylabel('TMB (mut/Mb)')
    ax.set_title('TMB by Immune Subtype')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    fig_path = "results/figures/tmb_by_subtype.png"
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"\n  ✓ Plot saved: {fig_path}")

# =============================================================================
# 6. TMB vs PD1/PD-L1 Signature
# =============================================================================

print("\n" + "=" * 60)
print(" TMB vs PD1/PD-L1 SIGNATURE")
print("=" * 60)

if 'PD1_PDL1_signature' in merged_df.columns:
    corr, p = spearmanr(merged_df['TMB_mut_per_Mb'], merged_df['PD1_PDL1_signature'])
    print(f"\n  Spearman correlation: r = {corr:.4f}, p = {p:.4f}")
    
    # Scatter plot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(merged_df['TMB_mut_per_Mb'], merged_df['PD1_PDL1_signature'], alpha=0.5)
    ax.set_xlabel('TMB (mut/Mb)')
    ax.set_ylabel('PD1/PD-L1 Signature')
    ax.set_title(f'TMB vs PD1/PD-L1 Signature\nr = {corr:.4f}, p = {p:.4f}')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    fig_path = "results/figures/tmb_vs_pd1.png"
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"  ✓ Plot saved: {fig_path}")

# =============================================================================
# 7. TMB vs IIR Group
# =============================================================================

print("\n" + "=" * 60)
print(" TMB vs IIR GROUP")
print("=" * 60)

if 'IIR_group' in merged_df.columns:
    print("\n  TMB by IIR group:")
    iir_data = {}
    for group in ['High_ICB_ready', 'Intermediate', 'Poor_ICB_ready']:
        subset = merged_df[merged_df['IIR_group'] == group]
        if len(subset) > 0:
            mean_val = subset['TMB_mut_per_Mb'].mean()
            std_val = subset['TMB_mut_per_Mb'].std()
            median_val = subset['TMB_mut_per_Mb'].median()
            print(f"    {group}: mean={mean_val:.2f} ± {std_val:.2f}, median={median_val:.2f}, n={len(subset)}")
            iir_data[group] = subset['TMB_mut_per_Mb'].values
    
    # Boxplot (fixed)
    fig, ax = plt.subplots(figsize=(8, 6))
    data = [iir_data.get(g, []) for g in ['High_ICB_ready', 'Intermediate', 'Poor_ICB_ready']]
    bp = ax.boxplot(data, patch_artist=True)
    ax.set_xticklabels(['High', 'Intermediate', 'Poor'])
    ax.set_ylabel('TMB (mut/Mb)')
    ax.set_title('TMB by IIR Group')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    fig_path = "results/figures/tmb_by_iir.png"
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"\n  ✓ Plot saved: {fig_path}")

# =============================================================================
# 8. TMB vs Survival
# =============================================================================

print("\n" + "=" * 60)
print(" TMB vs SURVIVAL")
print("=" * 60)

if 'os_time' in merged_df.columns and 'os_event' in merged_df.columns:
    # Clean survival data
    merged_df['os_time'] = pd.to_numeric(merged_df['os_time'], errors='coerce')
    merged_df['os_event'] = pd.to_numeric(merged_df['os_event'], errors='coerce')
    
    surv_df = merged_df.dropna(subset=['os_time', 'os_event', 'TMB_mut_per_Mb'])
    print(f"\n  Survival data: {len(surv_df)} samples")
    print(f"  Events: {surv_df['os_event'].sum()}")
    
    # TMB tertiles
    tertiles = np.percentile(surv_df['TMB_mut_per_Mb'], [33.33, 66.67])
    surv_df['TMB_tertile'] = pd.cut(
        surv_df['TMB_mut_per_Mb'],
        bins=[-np.inf, tertiles[0], tertiles[1], np.inf],
        labels=['Low', 'Mid', 'High']
    )
    print(f"\n  TMB tertiles: {tertiles[0]:.2f}, {tertiles[1]:.2f}")
    print(surv_df['TMB_tertile'].value_counts().sort_index())
    
    # KM Curves
    fig, ax = plt.subplots(figsize=(10, 6))
    kmf = KaplanMeierFitter()
    
    for group in ['Low', 'Mid', 'High']:
        mask = surv_df['TMB_tertile'] == group
        if mask.sum() > 0:
            kmf.fit(
                durations=surv_df.loc[mask, 'os_time'],
                event_observed=surv_df.loc[mask, 'os_event'],
                label=f'TMB {group}'
            )
            kmf.plot_survival_function(ax=ax, ci_show=True)
    
    ax.set_title('Survival by TMB Tertile (Standardized)', fontsize=14)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Survival Probability')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    fig_path = "results/figures/tmb_survival_km.png"
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"\n  ✓ KM curve saved: {fig_path}")
    
    # Log-rank test
    low_mask = surv_df['TMB_tertile'] == 'Low'
    high_mask = surv_df['TMB_tertile'] == 'High'
    
    if low_mask.sum() > 0 and high_mask.sum() > 0:
        logrank = logrank_test(
            surv_df.loc[low_mask, 'os_time'],
            surv_df.loc[high_mask, 'os_time'],
            surv_df.loc[low_mask, 'os_event'],
            surv_df.loc[high_mask, 'os_event']
        )
        print(f"\n  Log-rank test (Low vs High): p = {logrank.p_value:.4f}")
    
    # Cox PH
    cph = CoxPHFitter()
    cox_df = surv_df[['os_time', 'os_event', 'TMB_mut_per_Mb']].dropna()
    cph.fit(cox_df, 'os_time', 'os_event')
    
    hr = np.exp(cph.params_['TMB_mut_per_Mb'])
    hr_lower = np.exp(cph.confidence_intervals_.loc['TMB_mut_per_Mb', '95% lower-bound'])
    hr_upper = np.exp(cph.confidence_intervals_.loc['TMB_mut_per_Mb', '95% upper-bound'])
    p_val = cph.summary.loc['TMB_mut_per_Mb', 'p']
    
    print(f"\n  Cox PH (TMB continuous):")
    print(f"    HR = {hr:.4f} (95% CI: {hr_lower:.4f} - {hr_upper:.4f})")
    print(f"    p = {p_val:.4f}")

# =============================================================================
# 9. CORRELATION MATRIX
# =============================================================================

print("\n" + "=" * 60)
print(" CORRELATION MATRIX")
print("=" * 60)

# Select key columns
corr_cols = ['TMB_mut_per_Mb', 'ImmuneScore', 'PD1_PDL1_signature', 'Immune_Exclusion_index']
available_cols = [c for c in corr_cols if c in merged_df.columns]

if len(available_cols) >= 2:
    corr_matrix = merged_df[available_cols].corr(method='spearman')
    
    print("\n  Spearman correlation matrix:")
    print(corr_matrix.round(4))
    
    # Heatmap
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', center=0, fmt='.3f', 
                square=True, ax=ax)
    ax.set_title('Correlation Matrix (Spearman)')
    
    plt.tight_layout()
    fig_path = "results/figures/correlation_matrix.png"
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"\n  ✓ Heatmap saved: {fig_path}")

# =============================================================================
# 10. SAVE RESULTS
# =============================================================================

print("\nSaving results...")

# Save merged data
merged_df.to_csv("results/phase_vi/phase_vi_standardized_tmb.tsv", sep='\t', index=False)
print(f"  ✓ Merged data saved: results/phase_vi/phase_vi_standardized_tmb.tsv")

# Create summary
summary = f"""PHASE VI CORRELATIONS — STANDARDIZED TMB
========================================
TMB median: {merged_df['TMB_mut_per_Mb'].median():.2f} mut/Mb
TMB range: {merged_df['TMB_mut_per_Mb'].min():.2f} - {merged_df['TMB_mut_per_Mb'].max():.2f} mut/Mb
Samples: {len(merged_df)}

Key Correlations (Spearman):
  TMB vs ImmuneScore: r = {spearmanr(merged_df['TMB_mut_per_Mb'], merged_df['ImmuneScore'])[0]:.4f}, p = {spearmanr(merged_df['TMB_mut_per_Mb'], merged_df['ImmuneScore'])[1]:.4f}
  TMB vs PD1/PD-L1: r = {spearmanr(merged_df['TMB_mut_per_Mb'], merged_df['PD1_PDL1_signature'])[0]:.4f}, p = {spearmanr(merged_df['TMB_mut_per_Mb'], merged_df['PD1_PDL1_signature'])[1]:.4f}

TMB by Immune Subtype (median):
  BLIS: {merged_df[merged_df['immune_subtype'] == 'BLIS']['TMB_mut_per_Mb'].median():.2f}
  IM: {merged_df[merged_df['immune_subtype'] == 'IM']['TMB_mut_per_Mb'].median():.2f}
  BLIA: {merged_df[merged_df['immune_subtype'] == 'BLIA']['TMB_mut_per_Mb'].median():.2f}

TMB by IIR Group (median):
  High_ICB_ready: {merged_df[merged_df['IIR_group'] == 'High_ICB_ready']['TMB_mut_per_Mb'].median():.2f}
  Intermediate: {merged_df[merged_df['IIR_group'] == 'Intermediate']['TMB_mut_per_Mb'].median():.2f}
  Poor_ICB_ready: {merged_df[merged_df['IIR_group'] == 'Poor_ICB_ready']['TMB_mut_per_Mb'].median():.2f}
"""

with open("results/phase_vi/phase_vi_summary.txt", "w") as f:
    f.write(summary)

print(f"  ✓ Summary saved: results/phase_vi/phase_vi_summary.txt")

# =============================================================================
# 11. SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print(" SUMMARY — PHASE VI CORRELATIONS")
print("=" * 60)

print(f"\n  TMB median: {merged_df['TMB_mut_per_Mb'].median():.2f} mut/Mb")
print(f"  Samples: {len(merged_df)}")

if 'ImmuneScore' in merged_df.columns:
    corr, p = spearmanr(merged_df['TMB_mut_per_Mb'], merged_df['ImmuneScore'])
    print(f"\n  TMB vs ImmuneScore: r = {corr:.4f}, p = {p:.4f}")

if 'immune_subtype' in merged_df.columns:
    print("\n  TMB by Immune Subtype (median):")
    for st in ['BLIS', 'IM', 'BLIA']:
        subset = merged_df[merged_df['immune_subtype'] == st]
        if len(subset) > 0:
            print(f"    {st}: {subset['TMB_mut_per_Mb'].median():.2f} (n={len(subset)})")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ PHASE VI CORRELATIONS COMPLETE")