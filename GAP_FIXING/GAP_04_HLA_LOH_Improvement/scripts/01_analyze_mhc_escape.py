#!/usr/bin/env python3
"""
Script: 01_analyze_mhc_escape.py
Purpose: Analyze MHC-I escape as primary immune escape metric
Author: Bhaskararao Ch (Baashi)
GAP: 04 — HLA LOH Proxy Improvement
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# 1. SETUP
# =============================================================================

os.chdir("D:/Baashi/TNBC_project/Immune_Spatial_Immunotherapy/GAP_FIXING/GAP_04_HLA_LOH_Improvement")

print("=" * 60)
print(" MHC-I ESCAPE ANALYSIS")
print(" Primary Immune Escape Metric Validation")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

os.makedirs("results", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# =============================================================================
# 2. LOAD DATA
# =============================================================================

print("Loading data...")

# Load HLA escape data from Phase X
escape_file = "../../Phase_X_HLA_ImmuneEscape/results/hla/HLA_escape_TNBC_like_TCGA.tsv"

if not os.path.exists(escape_file):
    escape_file = "../../../Immune_Spatial_Immunotherapy/Phase_X_HLA_ImmuneEscape/results/hla/HLA_escape_TNBC_like_TCGA.tsv"

if not os.path.exists(escape_file):
    print(f"ERROR: Escape file not found: {escape_file}")
    exit(1)

df = pd.read_csv(escape_file, sep='\t')
print(f"  Loaded: {df.shape[0]} samples")
print(f"  Columns: {df.columns.tolist()}")

# =============================================================================
# 3. CHECK MHC-I ESCAPE GROUPS
# =============================================================================

print("\n" + "=" * 60)
print(" MHC-I ESCAPE GROUPS")
print("=" * 60)

# Check MHC-I tertile
if 'MHC_I_axis_tertile' in df.columns:
    print(f"\n  MHC-I tertile distribution:")
    print(df['MHC_I_axis_tertile'].value_counts().sort_index())

# Check MHC-I low flag
if 'MHC_I_low_flag' in df.columns:
    print(f"\n  MHC-I low flag (0=high, 1=low):")
    print(df['MHC_I_low_flag'].value_counts())

# Check HLA escape group
if 'HLA_escape_group' in df.columns:
    print(f"\n  HLA escape group:")
    print(df['HLA_escape_group'].value_counts())

# Check HLA LOH (for comparison - will drop)
if 'HLA_LOH_proxy_group' in df.columns:
    print(f"\n  HLA LOH proxy group (for reference - WILL DROP):")
    print(df['HLA_LOH_proxy_group'].value_counts())

# =============================================================================
# 4. SURVIVAL ANALYSIS: MHC-I ESCAPE GROUPS
# =============================================================================

print("\n" + "=" * 60)
print(" SURVIVAL ANALYSIS: MHC-I ESCAPE GROUPS")
print("=" * 60)

# Prepare survival data
if 'os_time' in df.columns and 'os_event' in df.columns:
    # Clean survival data
    df['os_time'] = pd.to_numeric(df['os_time'], errors='coerce')
    df['os_event'] = pd.to_numeric(df['os_event'], errors='coerce')
    
    # Remove NA
    surv_df = df.dropna(subset=['os_time', 'os_event', 'HLA_escape_group'])
    print(f"\n  Survival data: {surv_df.shape[0]} samples")
    print(f"  Events: {surv_df['os_event'].sum()}")
    
    # KM Curves
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: HLA escape group
    ax = axes[0]
    kmf = KaplanMeierFitter()
    
    for group in ['No_escape', 'Partial_escape']:
        mask = surv_df['HLA_escape_group'] == group
        if mask.sum() > 0:
            kmf.fit(
                durations=surv_df.loc[mask, 'os_time'],
                event_observed=surv_df.loc[mask, 'os_event'],
                label=f'{group}'
            )
            kmf.plot_survival_function(ax=ax, ci_show=True)
    
    ax.set_title('Survival by HLA Escape Group', fontsize=14)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Survival Probability')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Log-rank test
    no_escape_mask = surv_df['HLA_escape_group'] == 'No_escape'
    partial_escape_mask = surv_df['HLA_escape_group'] == 'Partial_escape'
    
    if no_escape_mask.sum() > 0 and partial_escape_mask.sum() > 0:
        logrank = logrank_test(
            surv_df.loc[no_escape_mask, 'os_time'],
            surv_df.loc[partial_escape_mask, 'os_time'],
            surv_df.loc[no_escape_mask, 'os_event'],
            surv_df.loc[partial_escape_mask, 'os_event']
        )
        print(f"\n  Log-rank test (No_escape vs Partial_escape):")
        print(f"    p-value: {logrank.p_value:.4f}")
        
        # Cox PH
        cox_df = surv_df[['os_time', 'os_event', 'HLA_escape_group']].copy()
        cox_df['is_partial_escape'] = (cox_df['HLA_escape_group'] == 'Partial_escape').astype(int)
        
        cph = CoxPHFitter()
        cph.fit(cox_df[['os_time', 'os_event', 'is_partial_escape']], 'os_time', 'os_event')
        
        hr = np.exp(cph.params_['is_partial_escape'])
        hr_lower = np.exp(cph.confidence_intervals_.loc['is_partial_escape', '95% lower-bound'])
        hr_upper = np.exp(cph.confidence_intervals_.loc['is_partial_escape', '95% upper-bound'])
        p_val = cph.summary.loc['is_partial_escape', 'p']
        
        print(f"\n  Cox PH (Partial_escape vs No_escape):")
        print(f"    HR = {hr:.4f} (95% CI: {hr_lower:.4f} - {hr_upper:.4f})")
        print(f"    p = {p_val:.4f}")
    
    # Plot 2: MHC-I tertile
    ax = axes[1]
    
    if 'MHC_I_axis_tertile' in surv_df.columns:
        for group in ['Low', 'Mid', 'High']:
            mask = surv_df['MHC_I_axis_tertile'] == group
            if mask.sum() > 0:
                kmf.fit(
                    durations=surv_df.loc[mask, 'os_time'],
                    event_observed=surv_df.loc[mask, 'os_event'],
                    label=f'MHC-I {group}'
                )
                kmf.plot_survival_function(ax=ax, ci_show=True)
        
        ax.set_title('Survival by MHC-I Expression', fontsize=14)
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('Survival Probability')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    fig_path = "results/figures/mhc_escape_survival.png"
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"\n  ✓ KM curves saved: {fig_path}")

# =============================================================================
# 5. MHC-I ESCAPE vs IIR GROUPS
# =============================================================================

print("\n" + "=" * 60)
print(" MHC-I ESCAPE vs IIR GROUPS")
print("=" * 60)

if 'IIR_group' in df.columns:
    print(f"\n  IIR group distribution:")
    print(df['IIR_group'].value_counts())
    
    # Cross-tabulation: IIR_group × HLA_escape_group
    crosstab = pd.crosstab(df['IIR_group'], df['HLA_escape_group'])
    print(f"\n  Cross-tabulation (IIR_group × HLA_escape_group):")
    print(crosstab)
    
    # Chi-square test
    chi2, p, dof, expected = chi2_contingency(crosstab)
    print(f"\n  Chi-square test:")
    print(f"    chi2 = {chi2:.4f}")
    print(f"    p = {p:.4f}")
    
    # Create 4-quadrant groups
    df['Escape_IIR_group'] = df['IIR_group'] + " + " + df['HLA_escape_group']
    print(f"\n  4-quadrant groups:")
    print(df['Escape_IIR_group'].value_counts())

# =============================================================================
# 6. SURVIVAL: 4-QUADRANT GROUPS
# =============================================================================

print("\n" + "=" * 60)
print(" SURVIVAL: 4-QUADRANT GROUPS")
print("=" * 60)

if 'Escape_IIR_group' in df.columns and 'os_time' in df.columns:
    surv_df = df.dropna(subset=['os_time', 'os_event', 'Escape_IIR_group'])
    
    # KM by 4-quadrant
    fig, ax = plt.subplots(figsize=(10, 6))
    kmf = KaplanMeierFitter()
    
    # Order groups for better visualization
    group_order = [
        'High_ICB_ready + No_escape',
        'High_ICB_ready + Partial_escape',
        'Poor_ICB_ready + No_escape',
        'Poor_ICB_ready + Partial_escape'
    ]
    
    for group in group_order:
        mask = surv_df['Escape_IIR_group'] == group
        if mask.sum() > 0:
            kmf.fit(
                durations=surv_df.loc[mask, 'os_time'],
                event_observed=surv_df.loc[mask, 'os_event'],
                label=group
            )
            kmf.plot_survival_function(ax=ax, ci_show=True)
    
    ax.set_title('Survival by 4-Quadrant Groups', fontsize=14)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Survival Probability')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    fig_path = "results/figures/4_quadrant_survival.png"
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"\n  ✓ 4-quadrant KM saved: {fig_path}")
    
    # Cox PH for 4-quadrant (High_ICB_ready + No_escape as reference)
    cox_df = surv_df[['os_time', 'os_event', 'Escape_IIR_group']].copy()
    
    # Create dummy variables
    for group in group_order:
        cox_df[f'is_{group.replace(" + ", "_")}'] = (cox_df['Escape_IIR_group'] == group).astype(int)
    
    # Reference: High_ICB_ready_No_escape
    ref = 'High_ICB_ready_No_escape'
    if ref in cox_df.columns:
        cox_df = cox_df.drop(columns=[ref])
        
        cph = CoxPHFitter()
        cph.fit(cox_df[['os_time', 'os_event'] + [c for c in cox_df.columns if c.startswith('is_')]], 
                'os_time', 'os_event')
        
        print(f"\n  Cox PH (Reference: {ref}):")
        for group in group_order:
            col = f'is_{group.replace(" + ", "_")}'
            if col in cph.params_.index:
                hr = np.exp(cph.params_[col])
                hr_lower = np.exp(cph.confidence_intervals_.loc[col, '95% lower-bound'])
                hr_upper = np.exp(cph.confidence_intervals_.loc[col, '95% upper-bound'])
                p_val = cph.summary.loc[col, 'p']
                print(f"    {group}: HR={hr:.4f} (95% CI: {hr_lower:.4f} - {hr_upper:.4f}), p={p_val:.4f}")

# =============================================================================
# 7. SAVE RESULTS
# =============================================================================

print("\nSaving results...")

# Save processed data
output_file = "results/mhc_escape_analysis_data.tsv"
df.to_csv(output_file, sep='\t', index=False)
print(f"  ✓ Data saved: {output_file}")

# Create summary
summary = f"""MHC-I ESCAPE ANALYSIS SUMMARY
========================================
Dataset: TCGA TNBC-like (n={len(df)})

MHC-I Escape Groups:
  No_escape: {(df['HLA_escape_group'] == 'No_escape').sum()}
  Partial_escape: {(df['HLA_escape_group'] == 'Partial_escape').sum()}

HLA LOH Proxy (DROPPED):
  Low: {(df['HLA_LOH_proxy_group'] == 'Low').sum()}
  High: {(df['HLA_LOH_proxy_group'] == 'High').sum()}
  NOTE: Only ~7 patients had HLA structural variants.

IIR Groups:
  High_ICB_ready: {(df['IIR_group'] == 'High_ICB_ready').sum()}
  Intermediate: {(df['IIR_group'] == 'Intermediate').sum()}
  Poor_ICB_ready: {(df['IIR_group'] == 'Poor_ICB_ready').sum()}

4-Quadrant Groups (IIR × Escape):
"""
for group in ['High_ICB_ready + No_escape', 'High_ICB_ready + Partial_escape',
              'Poor_ICB_ready + No_escape', 'Poor_ICB_ready + Partial_escape']:
    count = (df['Escape_IIR_group'] == group).sum()
    summary += f"  {group}: {count}\n"

with open("results/mhc_escape_summary.txt", "w") as f:
    f.write(summary)

print(f"  ✓ Summary saved: results/mhc_escape_summary.txt")

# =============================================================================
# 8. SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print(" SUMMARY — MHC-I ESCAPE ANALYSIS")
print("=" * 60)

print(f"\n  Total samples: {len(df)}")
print(f"\n  MHC-I Escape:")
print(f"    No_escape: {(df['HLA_escape_group'] == 'No_escape').sum()}")
print(f"    Partial_escape: {(df['HLA_escape_group'] == 'Partial_escape').sum()}")

print(f"\n  HLA LOH Proxy (DROPPED):")
print(f"    Low: {(df['HLA_LOH_proxy_group'] == 'Low').sum()}")
print(f"    High: {(df['HLA_LOH_proxy_group'] == 'High').sum()}")
print(f"    NOTE: Only ~7 patients had HLA structural variants.")

if 'Escape_IIR_group' in df.columns:
    print(f"\n  4-Quadrant Groups:")
    for group in ['High_ICB_ready + No_escape', 'High_ICB_ready + Partial_escape',
                  'Poor_ICB_ready + No_escape', 'Poor_ICB_ready + Partial_escape']:
        print(f"    {group}: {(df['Escape_IIR_group'] == group).sum()}")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ MHC-I ESCAPE ANALYSIS COMPLETE")