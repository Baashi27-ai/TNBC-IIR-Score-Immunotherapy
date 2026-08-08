#!/usr/bin/env python3
"""
Script: 07_validate_iir.py
Purpose: Validate IIR score with survival analysis on METABRIC
Author: Bhaskararao Ch (Baashi)
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# 1. SETUP
# =============================================================================

os.chdir("D:/Baashi/TNBC_project/Immune_Spatial_Immunotherapy/GAP_FIXING/GAP_01_External_Validation_METABRIC")

print("=" * 60)
print(" IIR SURVIVAL VALIDATION - METABRIC TNBC")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

# Create output directories
os.makedirs("results/survival", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# =============================================================================
# 2. LOAD INPUTS
# =============================================================================

print("Loading inputs...")

# Load IIR scores
iir_file = "results/iir/iir_score_metabric.tsv"
if not os.path.exists(iir_file):
    print(f"ERROR: IIR file not found: {iir_file}")
    print("Please run 06_compute_iir.py first.")
    exit(1)

iir_df = pd.read_csv(iir_file, sep='\t')
print(f"  IIR scores: {iir_df.shape[0]} samples")

# Load clinical data - FIXED PATH
clinical_file = "D:/Baashi/TNBC_project/M9_external_validation/inputs/METABRIC_TNBC_clinical_with_survival.tsv"
if not os.path.exists(clinical_file):
    print(f"ERROR: Clinical file not found: {clinical_file}")
    print("Trying alternative path...")
    clinical_file = "../../../M9_external_validation/inputs/METABRIC_TNBC_clinical_with_survival.tsv"
    if not os.path.exists(clinical_file):
        print(f"ERROR: Clinical file still not found: {clinical_file}")
        exit(1)

clinical_df = pd.read_csv(clinical_file, sep='\t')
print(f"  Clinical data: {clinical_df.shape[0]} samples")

# =============================================================================
# 3. MERGE DATA
# =============================================================================

print("\nMerging data...")

# Merge on sample_id
merged = pd.merge(iir_df, clinical_df, left_on='sample_id', right_on='SAMPLE_ID')

# Keep only necessary columns
merged = merged[[
    'sample_id',
    'IIR_score',
    'IIR_group',
    'OS_MONTHS',
    'OS_STATUS'
]].copy()

# Clean survival data
merged['OS_MONTHS'] = pd.to_numeric(merged['OS_MONTHS'], errors='coerce')
merged['OS_STATUS'] = merged['OS_STATUS'].apply(
    lambda x: 1 if 'DECEASED' in str(x) else 0
)

# Remove rows with missing survival
merged = merged.dropna(subset=['OS_MONTHS', 'OS_STATUS'])
print(f"  Merged: {merged.shape[0]} samples with survival data")

# =============================================================================
# 4. SURVIVAL ANALYSIS
# =============================================================================

print("\n" + "=" * 60)
print(" SURVIVAL ANALYSIS")
print("=" * 60)

# 4.1 Kaplan-Meier Curves
print("\nGenerating Kaplan-Meier curves...")

# Create figure
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: IIR groups (High vs Intermediate vs Poor)
ax = axes[0]
kmf = KaplanMeierFitter()

for group in ['High_ICB_ready', 'Intermediate', 'Poor_ICB_ready']:
    mask = merged['IIR_group'] == group
    if mask.sum() > 0:
        kmf.fit(
            durations=merged.loc[mask, 'OS_MONTHS'],
            event_observed=merged.loc[mask, 'OS_STATUS'],
            label=group
        )
        kmf.plot_survival_function(ax=ax, ci_show=True)

ax.set_title('IIR Groups - Overall Survival', fontsize=14)
ax.set_xlabel('Time (months)')
ax.set_ylabel('Survival Probability')
ax.legend(title='IIR Group')
ax.grid(True, alpha=0.3)

# Log-rank test for High vs Intermediate
mask_high = merged['IIR_group'] == 'High_ICB_ready'
mask_inter = merged['IIR_group'] == 'Intermediate'

if mask_high.sum() > 0 and mask_inter.sum() > 0:
    results = logrank_test(
        merged.loc[mask_high, 'OS_MONTHS'],
        merged.loc[mask_inter, 'OS_MONTHS'],
        merged.loc[mask_high, 'OS_STATUS'],
        merged.loc[mask_inter, 'OS_STATUS']
    )
    p_high_vs_inter = results.p_value
    print(f"  Log-rank p-value (High vs Intermediate): {p_high_vs_inter:.4f}")

# Log-rank test for High vs Poor
mask_poor = merged['IIR_group'] == 'Poor_ICB_ready'
if mask_high.sum() > 0 and mask_poor.sum() > 0:
    results = logrank_test(
        merged.loc[mask_high, 'OS_MONTHS'],
        merged.loc[mask_poor, 'OS_MONTHS'],
        merged.loc[mask_high, 'OS_STATUS'],
        merged.loc[mask_poor, 'OS_STATUS']
    )
    p_high_vs_poor = results.p_value
    print(f"  Log-rank p-value (High vs Poor): {p_high_vs_poor:.4f}")

# Plot 2: IIR tertiles (Low vs Mid vs High)
ax = axes[1]
# Recreate tertiles from IIR_score
tertiles = np.percentile(merged['IIR_score'], [33.33, 66.67])
merged['IIR_tertile'] = pd.cut(
    merged['IIR_score'],
    bins=[-np.inf, tertiles[0], tertiles[1], np.inf],
    labels=['Low', 'Mid', 'High']
)

for group in ['High', 'Mid', 'Low']:
    mask = merged['IIR_tertile'] == group
    if mask.sum() > 0:
        kmf.fit(
            durations=merged.loc[mask, 'OS_MONTHS'],
            event_observed=merged.loc[mask, 'OS_STATUS'],
            label=f'IIR {group}'
        )
        kmf.plot_survival_function(ax=ax, ci_show=True)

ax.set_title('IIR Tertiles - Overall Survival', fontsize=14)
ax.set_xlabel('Time (months)')
ax.set_ylabel('Survival Probability')
ax.legend(title='IIR Tertile')
ax.grid(True, alpha=0.3)

plt.tight_layout()
fig_path = "results/figures/iir_km_curve_metabric.png"
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
print(f"  ✓ KM curves saved to: {fig_path}")

# =============================================================================
# 5. COX PROPORTIONAL HAZARDS
# =============================================================================

print("\nRunning Cox Proportional Hazards models...")

cox_results = []

# Model 1: IIR_group (High vs Intermediate)
print("\n  Model 1: IIR_group (Reference: Intermediate)")
cox_df = merged[merged['IIR_group'].isin(['High_ICB_ready', 'Intermediate'])].copy()
if len(cox_df) > 0:
    cox_df['IIR_high'] = (cox_df['IIR_group'] == 'High_ICB_ready').astype(int)
    cph = CoxPHFitter()
    cph.fit(cox_df[['OS_MONTHS', 'OS_STATUS', 'IIR_high']], 'OS_MONTHS', 'OS_STATUS')
    
    # Extract results
    hr = np.exp(cph.params_['IIR_high'])
    ci_lower = np.exp(cph.confidence_intervals_.loc['IIR_high', '95% lower-bound'])
    ci_upper = np.exp(cph.confidence_intervals_.loc['IIR_high', '95% upper-bound'])
    p_val = cph.summary.loc['IIR_high', 'p']
    
    cox_results.append({
        'model': 'IIR_High_vs_Intermediate',
        'term': 'High_ICB_ready',
        'HR': hr,
        'HR_lower': ci_lower,
        'HR_upper': ci_upper,
        'p_value': p_val,
        'n': len(cox_df)
    })
    
    print(f"    HR = {hr:.4f} (95% CI: {ci_lower:.4f} - {ci_upper:.4f})")
    print(f"    p = {p_val:.4f}")

# Model 2: IIR_group (Poor vs Intermediate)
print("\n  Model 2: IIR_group (Reference: Intermediate)")
cox_df = merged[merged['IIR_group'].isin(['Poor_ICB_ready', 'Intermediate'])].copy()
if len(cox_df) > 0:
    cox_df['IIR_poor'] = (cox_df['IIR_group'] == 'Poor_ICB_ready').astype(int)
    cph = CoxPHFitter()
    cph.fit(cox_df[['OS_MONTHS', 'OS_STATUS', 'IIR_poor']], 'OS_MONTHS', 'OS_STATUS')
    
    hr = np.exp(cph.params_['IIR_poor'])
    ci_lower = np.exp(cph.confidence_intervals_.loc['IIR_poor', '95% lower-bound'])
    ci_upper = np.exp(cph.confidence_intervals_.loc['IIR_poor', '95% upper-bound'])
    p_val = cph.summary.loc['IIR_poor', 'p']
    
    cox_results.append({
        'model': 'IIR_Poor_vs_Intermediate',
        'term': 'Poor_ICB_ready',
        'HR': hr,
        'HR_lower': ci_lower,
        'HR_upper': ci_upper,
        'p_value': p_val,
        'n': len(cox_df)
    })
    
    print(f"    HR = {hr:.4f} (95% CI: {ci_lower:.4f} - {ci_upper:.4f})")
    print(f"    p = {p_val:.4f}")

# Model 3: IIR_score (continuous)
print("\n  Model 3: IIR_score (continuous)")
cph = CoxPHFitter()
cph.fit(merged[['OS_MONTHS', 'OS_STATUS', 'IIR_score']], 'OS_MONTHS', 'OS_STATUS')

hr = np.exp(cph.params_['IIR_score'])
ci_lower = np.exp(cph.confidence_intervals_.loc['IIR_score', '95% lower-bound'])
ci_upper = np.exp(cph.confidence_intervals_.loc['IIR_score', '95% upper-bound'])
p_val = cph.summary.loc['IIR_score', 'p']

cox_results.append({
    'model': 'IIR_score_continuous',
    'term': 'IIR_score',
    'HR': hr,
    'HR_lower': ci_lower,
    'HR_upper': ci_upper,
    'p_value': p_val,
    'n': len(merged)
})

print(f"    HR = {hr:.4f} (95% CI: {ci_lower:.4f} - {ci_upper:.4f})")
print(f"    p = {p_val:.4f}")

# =============================================================================
# 6. SAVE COX RESULTS
# =============================================================================

cox_df = pd.DataFrame(cox_results)
output_file = "results/survival/iir_survival_results.tsv"
cox_df.to_csv(output_file, sep='\t', index=False)
print(f"\n  ✓ Cox results saved to: {output_file}")

# =============================================================================
# 7. SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print(" SUMMARY")
print("=" * 60)

print(f"\nSurvival data summary:")
print(f"  Total samples: {len(merged)}")
print(f"  Events (deaths): {merged['OS_STATUS'].sum()}")
print(f"  Median follow-up: {merged['OS_MONTHS'].median():.1f} months")

print(f"\nIIR group survival:")
for group in ['High_ICB_ready', 'Intermediate', 'Poor_ICB_ready']:
    mask = merged['IIR_group'] == group
    if mask.sum() > 0:
        events = merged.loc[mask, 'OS_STATUS'].sum()
        median_os = merged.loc[mask, 'OS_MONTHS'].median()
        print(f"  {group}: n={mask.sum()}, events={events}, median OS={median_os:.1f} months")

print(f"\nCox model results (key):")
for row in cox_results:
    if row['model'] == 'IIR_High_vs_Intermediate':
        print(f"  High vs Intermediate: HR={row['HR']:.4f}, p={row['p_value']:.4f}")
    elif row['model'] == 'IIR_Poor_vs_Intermediate':
        print(f"  Poor vs Intermediate: HR={row['HR']:.4f}, p={row['p_value']:.4f}")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ IIR SURVIVAL VALIDATION COMPLETE")
print("   Next: 08_validate_ushape.py")