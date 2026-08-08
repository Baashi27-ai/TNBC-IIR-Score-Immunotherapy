#!/usr/bin/env python3
"""
Script: 01_compare_biomarkers.py
Purpose: Compare IIR score to existing clinical biomarkers
Author: Bhaskararao Ch (Baashi)
GAP: 06 — Biomarker Comparison
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# 1. SETUP
# =============================================================================

os.chdir("D:/Baashi/TNBC_project/Immune_Spatial_Immunotherapy/GAP_FIXING/GAP_06_Biomarker_Comparison")

print("=" * 60)
print(" BIOMARKER COMPARISON")
print(" IIR vs PD-L1 vs ImmuneScore vs TMB")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

os.makedirs("results", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# =============================================================================
# 2. LOAD DATA
# =============================================================================

print("Loading data...")

# Load IIR scores from GAP_01
iir_file = "../GAP_01_External_Validation_METABRIC/results/iir/iir_score_metabric.tsv"

if not os.path.exists(iir_file):
    iir_file = "../../GAP_01_External_Validation_METABRIC/results/iir/iir_score_metabric.tsv"

if not os.path.exists(iir_file):
    print(f"ERROR: IIR file not found: {iir_file}")
    exit(1)

iir_df = pd.read_csv(iir_file, sep='\t')
print(f"  IIR scores: {len(iir_df)} samples")
print(f"  IIR columns: {iir_df.columns.tolist()}")

# Load PD1/PD-L1 signature from GAP_01
pd1_file = "../GAP_01_External_Validation_METABRIC/results/immune_scores/pd1_signature_metabric.tsv"

if not os.path.exists(pd1_file):
    pd1_file = "../../GAP_01_External_Validation_METABRIC/results/immune_scores/pd1_signature_metabric.tsv"

if not os.path.exists(pd1_file):
    print(f"ERROR: PD1 file not found: {pd1_file}")
    exit(1)

pd1_df = pd.read_csv(pd1_file, sep='\t')
print(f"  PD1/PD-L1 signature: {len(pd1_df)} samples")
print(f"  PD1 columns: {pd1_df.columns.tolist()}")

# Load ImmuneScore from GAP_01
immune_file = "../GAP_01_External_Validation_METABRIC/results/immune_scores/immune_scores_metabric.tsv"

if not os.path.exists(immune_file):
    immune_file = "../../GAP_01_External_Validation_METABRIC/results/immune_scores/immune_scores_metabric.tsv"

if not os.path.exists(immune_file):
    print(f"ERROR: ImmuneScore file not found: {immune_file}")
    exit(1)

immune_df = pd.read_csv(immune_file, sep='\t')
print(f"  ImmuneScore: {len(immune_df)} samples")
print(f"  Immune columns: {immune_df.columns.tolist()}")

# Load clinical data
clinical_file = "D:/Baashi/TNBC_project/M9_external_validation/inputs/METABRIC_TNBC_clinical_with_survival.tsv"

if not os.path.exists(clinical_file):
    print(f"ERROR: Clinical file not found: {clinical_file}")
    exit(1)

clinical_df = pd.read_csv(clinical_file, sep='\t')
print(f"  Clinical data: {len(clinical_df)} samples")

# =============================================================================
# 3. MERGE DATA
# =============================================================================

print("\nMerging data...")

# Start with IIR (columns: sample_id, IIR_score, IIR_group)
merged_df = iir_df[['sample_id', 'IIR_score', 'IIR_group']].copy()
merged_df.rename(columns={'IIR_score': 'IIR_score_norm'}, inplace=True)

# Add PD1 signature (columns: sample_id, PD1_PDL1_signature, PD1_group)
pd1_cols = ['sample_id', 'PD1_PDL1_signature']
if 'PD1_group' in pd1_df.columns:
    pd1_cols.append('PD1_group')
merged_df = pd.merge(merged_df, pd1_df[pd1_cols], on='sample_id', how='left')

# Add ImmuneScore (columns: sample_id, ImmuneScore_norm, immune_group)
immune_cols = ['sample_id', 'ImmuneScore_norm']
if 'immune_group' in immune_df.columns:
    immune_cols.append('immune_group')
merged_df = pd.merge(merged_df, immune_df[immune_cols], on='sample_id', how='left')

# Add clinical data
merged_df = pd.merge(merged_df, clinical_df[['SAMPLE_ID', 'OS_MONTHS', 'OS_STATUS', 'TMB_NONSYNONYMOUS']], 
                     left_on='sample_id', right_on='SAMPLE_ID', how='inner')

# Clean survival data
merged_df['OS_MONTHS'] = pd.to_numeric(merged_df['OS_MONTHS'], errors='coerce')
merged_df['OS_STATUS'] = merged_df['OS_STATUS'].apply(
    lambda x: 1 if 'DECEASED' in str(x) else 0
)

# Clean TMB
merged_df['TMB'] = pd.to_numeric(merged_df['TMB_NONSYNONYMOUS'], errors='coerce')

# Drop NA
merged_df = merged_df.dropna(subset=['OS_MONTHS', 'OS_STATUS'])
print(f"  Merged: {len(merged_df)} samples")

# =============================================================================
# 4. CORRELATION BETWEEN BIOMARKERS
# =============================================================================

print("\n" + "=" * 60)
print(" CORRELATION BETWEEN BIOMARKERS")
print("=" * 60)

biomarkers = ['IIR_score_norm', 'PD1_PDL1_signature', 'ImmuneScore_norm', 'TMB']
available_biomarkers = [b for b in biomarkers if b in merged_df.columns]
corr_matrix = merged_df[available_biomarkers].corr(method='spearman')

print("\n  Spearman correlation matrix:")
print(corr_matrix.round(4))

if len(available_biomarkers) >= 2:
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', center=0, fmt='.3f', square=True, ax=ax)
    ax.set_title('Biomarker Correlations (Spearman)')
    plt.tight_layout()
    fig_path = "results/figures/biomarker_correlations.png"
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"\n  ✓ Heatmap saved: {fig_path}")

# =============================================================================
# 5. COX MODELS
# =============================================================================

print("\n" + "=" * 60)
print(" COX PROPORTIONAL HAZARDS MODELS")
print("=" * 60)

# Prepare data for Cox
cox_vars = ['IIR_score_norm', 'PD1_PDL1_signature', 'ImmuneScore_norm', 'TMB']
available_vars = [v for v in cox_vars if v in merged_df.columns]

cox_data = merged_df[['OS_MONTHS', 'OS_STATUS'] + available_vars].dropna()

print(f"\n  Cox data: {len(cox_data)} samples")
print(f"  Events: {cox_data['OS_STATUS'].sum()}")

# Store results
model_results = []

# ---- Model 1: PD-L1 alone ----
if 'PD1_PDL1_signature' in available_vars:
    print("\n  Model 1: PD-L1 alone")
    cph1 = CoxPHFitter()
    cph1.fit(cox_data[['OS_MONTHS', 'OS_STATUS', 'PD1_PDL1_signature']], 'OS_MONTHS', 'OS_STATUS')
    hr1 = np.exp(cph1.params_['PD1_PDL1_signature'])
    p1 = cph1.summary.loc['PD1_PDL1_signature', 'p']
    print(f"    HR: {hr1:.4f}, p: {p1:.4f}")
    print(f"    AIC: {cph1.AIC_partial_:.2f}")
    print(f"    C-index: {cph1.concordance_index_:.4f}")
    model_results.append({
        'Model': 'PD-L1 alone',
        'AIC': cph1.AIC_partial_,
        'C-index': cph1.concordance_index_,
        'LogLik': cph1.log_likelihood_,
        'HR': hr1,
        'p': p1
    })

# ---- Model 2: ImmuneScore alone ----
if 'ImmuneScore_norm' in available_vars:
    print("\n  Model 2: ImmuneScore alone")
    cph2 = CoxPHFitter()
    cph2.fit(cox_data[['OS_MONTHS', 'OS_STATUS', 'ImmuneScore_norm']], 'OS_MONTHS', 'OS_STATUS')
    hr2 = np.exp(cph2.params_['ImmuneScore_norm'])
    p2 = cph2.summary.loc['ImmuneScore_norm', 'p']
    print(f"    HR: {hr2:.4f}, p: {p2:.4f}")
    print(f"    AIC: {cph2.AIC_partial_:.2f}")
    print(f"    C-index: {cph2.concordance_index_:.4f}")
    model_results.append({
        'Model': 'ImmuneScore alone',
        'AIC': cph2.AIC_partial_,
        'C-index': cph2.concordance_index_,
        'LogLik': cph2.log_likelihood_,
        'HR': hr2,
        'p': p2
    })

# ---- Model 3: TMB alone ----
if 'TMB' in available_vars:
    print("\n  Model 3: TMB alone")
    cox_data_tmb = cox_data.dropna(subset=['TMB'])
    if len(cox_data_tmb) > 10:
        cph3 = CoxPHFitter()
        cph3.fit(cox_data_tmb[['OS_MONTHS', 'OS_STATUS', 'TMB']], 'OS_MONTHS', 'OS_STATUS')
        hr3 = np.exp(cph3.params_['TMB'])
        p3 = cph3.summary.loc['TMB', 'p']
        print(f"    HR: {hr3:.4f}, p: {p3:.4f}")
        print(f"    AIC: {cph3.AIC_partial_:.2f}")
        print(f"    C-index: {cph3.concordance_index_:.4f}")
        model_results.append({
            'Model': 'TMB alone',
            'AIC': cph3.AIC_partial_,
            'C-index': cph3.concordance_index_,
            'LogLik': cph3.log_likelihood_,
            'HR': hr3,
            'p': p3
        })
    else:
        print("  ⚠️ TMB data insufficient for Cox model")

# ---- Model 4: IIR score ----
if 'IIR_score_norm' in available_vars:
    print("\n  Model 4: IIR score (full model)")
    cph4 = CoxPHFitter()
    cph4.fit(cox_data[['OS_MONTHS', 'OS_STATUS', 'IIR_score_norm']], 'OS_MONTHS', 'OS_STATUS')
    hr4 = np.exp(cph4.params_['IIR_score_norm'])
    p4 = cph4.summary.loc['IIR_score_norm', 'p']
    print(f"    HR: {hr4:.4f}, p: {p4:.4f}")
    print(f"    AIC: {cph4.AIC_partial_:.2f}")
    print(f"    C-index: {cph4.concordance_index_:.4f}")
    model_results.append({
        'Model': 'IIR alone',
        'AIC': cph4.AIC_partial_,
        'C-index': cph4.concordance_index_,
        'LogLik': cph4.log_likelihood_,
        'HR': hr4,
        'p': p4
    })

# ---- Model 5: IIR + PD-L1 ----
if 'IIR_score_norm' in available_vars and 'PD1_PDL1_signature' in available_vars:
    print("\n  Model 5: IIR + PD-L1 (combined)")
    cph5 = CoxPHFitter()
    cph5.fit(cox_data[['OS_MONTHS', 'OS_STATUS', 'IIR_score_norm', 'PD1_PDL1_signature']], 'OS_MONTHS', 'OS_STATUS')
    print(f"    AIC: {cph5.AIC_partial_:.2f}")
    print(f"    C-index: {cph5.concordance_index_:.4f}")
    print(f"    Log-likelihood: {cph5.log_likelihood_:.2f}")
    for var in ['IIR_score_norm', 'PD1_PDL1_signature']:
        hr = np.exp(cph5.params_[var])
        p = cph5.summary.loc[var, 'p']
        print(f"    {var}: HR={hr:.4f}, p={p:.4f}")
    model_results.append({
        'Model': 'IIR + PD-L1',
        'AIC': cph5.AIC_partial_,
        'C-index': cph5.concordance_index_,
        'LogLik': cph5.log_likelihood_,
        'HR': None,
        'p': None
    })

# =============================================================================
# 6. COMPARE MODELS
# =============================================================================

print("\n" + "=" * 60)
print(" MODEL COMPARISON")
print("=" * 60)

results_df = pd.DataFrame(model_results)
if len(results_df) > 0:
    print("\n  Model comparison:")
    print(results_df.round(4))

# =============================================================================
# 7. VISUALIZE MODEL COMPARISON
# =============================================================================

if len(results_df) >= 2:
    print("\n" + "=" * 60)
    print(" VISUALIZING MODEL COMPARISON")
    print("=" * 60)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # AIC
    ax = axes[0]
    colors = ['lightgray' if 'IIR' not in m else 'lightcoral' for m in results_df['Model']]
    bars = ax.bar(results_df['Model'], results_df['AIC'], color=colors)
    ax.set_ylabel('AIC (lower is better)')
    ax.set_title('AIC Comparison')
    ax.tick_params(axis='x', rotation=45)

    # C-index
    ax = axes[1]
    colors = ['lightgray' if 'IIR' not in m else 'lightcoral' for m in results_df['Model']]
    bars = ax.bar(results_df['Model'], results_df['C-index'], color=colors)
    ax.axhline(y=0.65, color='red', linestyle='--', label='Target C-index > 0.65')
    ax.set_ylabel('C-index (higher is better)')
    ax.set_title('C-index Comparison')
    ax.tick_params(axis='x', rotation=45)
    ax.legend()

    plt.tight_layout()
    fig_path = "results/figures/model_comparison.png"
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"  ✓ Model comparison plot saved: {fig_path}")

# =============================================================================
# 8. SAVE RESULTS
# =============================================================================

print("\nSaving results...")

# Save results
results_df.to_csv("results/biomarker_comparison_results.tsv", sep='\t', index=False)
print(f"  ✓ Results saved: results/biomarker_comparison_results.tsv")

# Save merged data
merged_df.to_csv("results/biomarker_data.tsv", sep='\t', index=False)
print(f"  ✓ Data saved: results/biomarker_data.tsv")

# =============================================================================
# 9. SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print(" SUMMARY — BIOMARKER COMPARISON")
print("=" * 60)

print(f"\n  Samples: {len(cox_data)}")
print(f"  Events: {cox_data['OS_STATUS'].sum()}")

print("\n  Model Comparison:")
for result in model_results:
    print(f"    {result['Model']}: AIC={result['AIC']:.2f}, C-index={result['C-index']:.4f}")

# Check success criteria
if 'IIR alone' in results_df['Model'].values:
    iir_row = results_df[results_df['Model'] == 'IIR alone'].iloc[0]
    print("\n  Success Criteria:")
    if iir_row['C-index'] > 0.65:
        print(f"    ✅ IIR C-index = {iir_row['C-index']:.4f} > 0.65")
    else:
        print(f"    ❌ IIR C-index = {iir_row['C-index']:.4f} (target > 0.65)")

    # Compare to PD-L1
    if 'PD-L1 alone' in results_df['Model'].values:
        pd1_row = results_df[results_df['Model'] == 'PD-L1 alone'].iloc[0]
        if iir_row['AIC'] < pd1_row['AIC']:
            print(f"    ✅ IIR AIC lower than PD-L1 (ΔAIC = {pd1_row['AIC'] - iir_row['AIC']:.2f})")
        else:
            print(f"    ❌ IIR AIC higher than PD-L1")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ BIOMARKER COMPARISON COMPLETE")