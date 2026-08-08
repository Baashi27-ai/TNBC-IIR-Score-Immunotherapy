#!/usr/bin/env python3
"""
Script: 08_validate_ushape.py
Purpose: Validate U-shaped survival pattern on METABRIC
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
print(" U-SHAPE SURVIVAL VALIDATION - METABRIC TNBC")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

os.makedirs("results/survival", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# =============================================================================
# 2. LOAD INPUTS
# =============================================================================

print("Loading inputs...")

# Load spatial metrics
spatial_file = "results/spatial_metrics/spatial_metrics_metabric.tsv"
spatial_df = pd.read_csv(spatial_file, sep='\t')
print(f"  Spatial metrics: {spatial_df.shape[0]} samples")

# Load clinical data
clinical_file = "D:/Baashi/TNBC_project/M9_external_validation/inputs/METABRIC_TNBC_clinical_with_survival.tsv"
clinical_df = pd.read_csv(clinical_file, sep='\t')
print(f"  Clinical data: {clinical_df.shape[0]} samples")

# =============================================================================
# 3. MERGE DATA
# =============================================================================

merged = pd.merge(spatial_df, clinical_df, left_on='sample_id', right_on='SAMPLE_ID')

# Clean survival
merged['OS_MONTHS'] = pd.to_numeric(merged['OS_MONTHS'], errors='coerce')
merged['OS_STATUS'] = merged['OS_STATUS'].apply(lambda x: 1 if 'DECEASED' in str(x) else 0)
merged = merged.dropna(subset=['OS_MONTHS', 'OS_STATUS'])

print(f"  Merged: {merged.shape[0]} samples")

# =============================================================================
# 4. U-SHAPE SURVIVAL ANALYSIS
# =============================================================================

print("\n" + "=" * 60)
print(" U-SHAPE SURVIVAL ANALYSIS")
print("=" * 60)

# Use Immune_Stroma_ratio tertiles
ratio_tertile = merged['Immune_Stroma_ratio_tertile']

# KM Curves
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: Ratio tertiles
ax = axes[0]
kmf = KaplanMeierFitter()

for group in ['Low', 'Mid', 'High']:
    mask = ratio_tertile == group
    if mask.sum() > 0:
        kmf.fit(
            durations=merged.loc[mask, 'OS_MONTHS'],
            event_observed=merged.loc[mask, 'OS_STATUS'],
            label=f'Ratio {group}'
        )
        kmf.plot_survival_function(ax=ax, ci_show=True)

ax.set_title('Immune-Stroma Ratio - Overall Survival', fontsize=14)
ax.set_xlabel('Time (months)')
ax.set_ylabel('Survival Probability')
ax.legend(title='Ratio Group')
ax.grid(True, alpha=0.3)

# Log-rank tests
mask_low = ratio_tertile == 'Low'
mask_mid = ratio_tertile == 'Mid'
mask_high = ratio_tertile == 'High'

if mask_low.sum() > 0 and mask_mid.sum() > 0:
    p_low_vs_mid = logrank_test(
        merged.loc[mask_low, 'OS_MONTHS'],
        merged.loc[mask_mid, 'OS_MONTHS'],
        merged.loc[mask_low, 'OS_STATUS'],
        merged.loc[mask_mid, 'OS_STATUS']
    ).p_value
    print(f"  Log-rank (Low vs Mid): p = {p_low_vs_mid:.4f}")

if mask_high.sum() > 0 and mask_mid.sum() > 0:
    p_high_vs_mid = logrank_test(
        merged.loc[mask_high, 'OS_MONTHS'],
        merged.loc[mask_mid, 'OS_MONTHS'],
        merged.loc[mask_high, 'OS_STATUS'],
        merged.loc[mask_mid, 'OS_STATUS']
    ).p_value
    print(f"  Log-rank (High vs Mid): p = {p_high_vs_mid:.4f}")

# Plot 2: Exclusion tertiles (inverse)
ax = axes[1]
excl_tertile = merged['Immune_Exclusion_index_tertile']

for group in ['Low', 'Mid', 'High']:
    mask = excl_tertile == group
    if mask.sum() > 0:
        kmf.fit(
            durations=merged.loc[mask, 'OS_MONTHS'],
            event_observed=merged.loc[mask, 'OS_STATUS'],
            label=f'Exclusion {group}'
        )
        kmf.plot_survival_function(ax=ax, ci_show=True)

ax.set_title('Immune Exclusion Index - Overall Survival', fontsize=14)
ax.set_xlabel('Time (months)')
ax.set_ylabel('Survival Probability')
ax.legend(title='Exclusion Group')
ax.grid(True, alpha=0.3)

plt.tight_layout()
fig_path = "results/figures/ushape_km_curve_metabric.png"
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
print(f"\n  ✓ KM curves saved to: {fig_path}")

# =============================================================================
# 5. COX MODELS
# =============================================================================

print("\nRunning Cox models...")

cox_results = []

# Model: Low vs Mid (Reference: Mid)
cox_df = merged[ratio_tertile.isin(['Low', 'Mid'])].copy()
if len(cox_df) > 0:
    cox_df['is_low'] = (cox_df['Immune_Stroma_ratio_tertile'] == 'Low').astype(int)
    cph = CoxPHFitter()
    cph.fit(cox_df[['OS_MONTHS', 'OS_STATUS', 'is_low']], 'OS_MONTHS', 'OS_STATUS')
    
    hr = np.exp(cph.params_['is_low'])
    ci_lower = np.exp(cph.confidence_intervals_.loc['is_low', '95% lower-bound'])
    ci_upper = np.exp(cph.confidence_intervals_.loc['is_low', '95% upper-bound'])
    p_val = cph.summary.loc['is_low', 'p']
    
    cox_results.append({
        'model': 'Low_vs_Mid',
        'HR': hr,
        'HR_lower': ci_lower,
        'HR_upper': ci_upper,
        'p_value': p_val,
        'n': len(cox_df)
    })
    print(f"  Low vs Mid: HR={hr:.4f}, p={p_val:.4f}")

# Model: High vs Mid (Reference: Mid)
cox_df = merged[ratio_tertile.isin(['High', 'Mid'])].copy()
if len(cox_df) > 0:
    cox_df['is_high'] = (cox_df['Immune_Stroma_ratio_tertile'] == 'High').astype(int)
    cph = CoxPHFitter()
    cph.fit(cox_df[['OS_MONTHS', 'OS_STATUS', 'is_high']], 'OS_MONTHS', 'OS_STATUS')
    
    hr = np.exp(cph.params_['is_high'])
    ci_lower = np.exp(cph.confidence_intervals_.loc['is_high', '95% lower-bound'])
    ci_upper = np.exp(cph.confidence_intervals_.loc['is_high', '95% upper-bound'])
    p_val = cph.summary.loc['is_high', 'p']
    
    cox_results.append({
        'model': 'High_vs_Mid',
        'HR': hr,
        'HR_lower': ci_lower,
        'HR_upper': ci_upper,
        'p_value': p_val,
        'n': len(cox_df)
    })
    print(f"  High vs Mid: HR={hr:.4f}, p={p_val:.4f}")

# =============================================================================
# 6. SAVE RESULTS
# =============================================================================

cox_df = pd.DataFrame(cox_results)
output_file = "results/survival/ushape_validation_results.tsv"
cox_df.to_csv(output_file, sep='\t', index=False)
print(f"\n  ✓ Cox results saved to: {output_file}")

# =============================================================================
# 7. SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print(" SUMMARY")
print("=" * 60)

print(f"\nU-shape validation summary:")
print(f"  Total samples: {len(merged)}")
print(f"  Events: {merged['OS_STATUS'].sum()}")

print(f"\nSurvival by Ratio group:")
for group in ['Low', 'Mid', 'High']:
    mask = ratio_tertile == group
    if mask.sum() > 0:
        events = merged.loc[mask, 'OS_STATUS'].sum()
        median_os = merged.loc[mask, 'OS_MONTHS'].median()
        print(f"  {group}: n={mask.sum()}, events={events}, median OS={median_os:.1f} months")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ U-SHAPE VALIDATION COMPLETE")