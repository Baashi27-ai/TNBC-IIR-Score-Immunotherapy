#!/usr/bin/env python3
"""
Script: 03_validate_iir_gse91061.py
Purpose: Validate IIR score on GSE91061 melanoma cohort (PD-1/CTLA-4 treated)
Author: Bhaskararao Ch (Baashi)
GAP: 02 — ICB-Treated Cohort Validation
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import zscore, fisher_exact
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# 1. SETUP
# =============================================================================

os.chdir("D:/Baashi/TNBC_project/Immune_Spatial_Immunotherapy/GAP_FIXING/GAP_02_ICB_Treated_Cohort_IMpassion")

print("=" * 60)
print(" GSE91061 VALIDATION — IIR SCORE")
print(" PD-1/CTLA-4 Treated Melanoma")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

os.makedirs("results/gse91061", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# =============================================================================
# 2. LOAD DATA
# =============================================================================

print("Loading GSE91061 data...")

# Load mapped expression
expr_file = "D:/Baashi/TNBC_project/ICB_Data/GSE91061/GSE91061_expression_mapped.csv"

if not os.path.exists(expr_file):
    print(f"ERROR: Expression file not found: {expr_file}")
    exit(1)

expr_df = pd.read_csv(expr_file)
print(f"  Expression: {expr_df.shape[0]} genes × {expr_df.shape[1]-1} samples")

# Load clinical data
clinical_file = "D:/Baashi/TNBC_project/ICB_Data/GSE91061/GSE91061_clean_clinical.csv"

if not os.path.exists(clinical_file):
    clinical_file = "D:/Baashi/TNBC_project/ICB_Data/GSE91061/GSE91061_clinical.csv"

if not os.path.exists(clinical_file):
    print(f"ERROR: Clinical file not found")
    exit(1)

clinical_df = pd.read_csv(clinical_file)
print(f"  Clinical: {clinical_df.shape[0]} samples")

# =============================================================================
# 3. PREPARE DATA
# =============================================================================

print("\nPreparing data...")

# Set gene symbols as index
expr_df = expr_df.set_index('Gene_Symbol')

# Get sample IDs (column names)
sample_ids = expr_df.columns.tolist()

# Transpose to samples × genes
expr_matrix = expr_df.T

print(f"  Expression matrix: {expr_matrix.shape[0]} samples × {expr_matrix.shape[1]} genes")

# =============================================================================
# 4. COMPUTE IIR SCORE
# =============================================================================

print("\nComputing IIR score...")

# PD1 signature genes
pd1_genes = [
    'PDCD1', 'CD274', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT', 'BTLA',
    'CXCL9', 'CXCL10', 'CXCL11', 'CXCL13',
    'IFNG', 'STAT1', 'IRF1',
    'CD8A', 'CD8B', 'GZMA', 'GZMB', 'PRF1', 'GNLY', 'NKG7',
    'CD3D', 'CD3E', 'CD3G', 'CD4',
    'HLA-DRA', 'HLA-DRB1', 'HLA-DQB1'
]

# Find available genes
available_genes = [g for g in pd1_genes if g in expr_matrix.columns]
print(f"  Available genes: {len(available_genes)}/{len(pd1_genes)}")

# Extract expression for available genes
gene_expr = expr_matrix[available_genes]

# Z-score normalize each gene
gene_z = gene_expr.apply(zscore, axis=0)

# PD1 signature (mean of z-scores)
expr_matrix['PD1_signature'] = gene_z.mean(axis=1)

# ImmuneScore (cytotoxic signature)
cytotoxic_genes = ['CD8A', 'CD8B', 'GZMA', 'GZMB', 'PRF1']
avail_cytotoxic = [g for g in cytotoxic_genes if g in available_genes]

if avail_cytotoxic:
    expr_matrix['ImmuneScore'] = expr_matrix[avail_cytotoxic].mean(axis=1)
    expr_matrix['ImmuneScore_norm'] = (expr_matrix['ImmuneScore'] - expr_matrix['ImmuneScore'].min()) / \
                                      (expr_matrix['ImmuneScore'].max() - expr_matrix['ImmuneScore'].min() + 1e-10)
else:
    expr_matrix['ImmuneScore_norm'] = expr_matrix['PD1_signature']

# DPI (simplified for GSE91061)
expr_matrix['DPI_norm'] = expr_matrix['ImmuneScore_norm']

# IIR Score (simplified: PD1 + ImmuneScore)
expr_matrix['IIR_score'] = (expr_matrix['PD1_signature'] + expr_matrix['ImmuneScore_norm']) / 2
expr_matrix['IIR_score_norm'] = (expr_matrix['IIR_score'] - expr_matrix['IIR_score'].min()) / \
                                (expr_matrix['IIR_score'].max() - expr_matrix['IIR_score'].min() + 1e-10)

# IIR group (median split)
median_iir = expr_matrix['IIR_score_norm'].median()
expr_matrix['IIR_group'] = np.where(expr_matrix['IIR_score_norm'] > median_iir, 'High', 'Low')

print(f"  IIR_score range: {expr_matrix['IIR_score_norm'].min():.4f} - {expr_matrix['IIR_score_norm'].max():.4f}")
print(f"  IIR median: {median_iir:.4f}")
print(f"  High group: {(expr_matrix['IIR_group'] == 'High').sum()}")
print(f"  Low group: {(expr_matrix['IIR_group'] == 'Low').sum()}")

# =============================================================================
# 5. MERGE WITH CLINICAL DATA
# =============================================================================

print("\nMerging with clinical data...")

# The clinical file has sample IDs
# We need to match them with expression sample IDs

# Check clinical columns
print(f"  Clinical columns: {clinical_df.columns.tolist()}")

# Try to find the right columns
if 'sample_id' in clinical_df.columns:
    sample_col = 'sample_id'
elif 'Sample' in clinical_df.columns:
    sample_col = 'Sample'
elif 'ID' in clinical_df.columns:
    sample_col = 'ID'
else:
    sample_col = clinical_df.columns[0]

print(f"  Using sample column: {sample_col}")

# Try to find response column
response_col = None
for col in clinical_df.columns:
    if 'response' in col.lower() or 'respond' in col.lower():
        response_col = col
        break

print(f"  Using response column: {response_col}")

# Merge
clinical_df = clinical_df.set_index(sample_col)
expr_matrix = expr_matrix.join(clinical_df, how='inner')

print(f"  Merged: {expr_matrix.shape[0]} samples")

# =============================================================================
# 6. VALIDATION: RESPONSE (ORR)
# =============================================================================

print("\n" + "=" * 60)
print(" VALIDATION: OBJECTIVE RESPONSE RATE (ORR)")
print("=" * 60)

response_results = {}

if response_col and response_col in expr_matrix.columns:
    # Map response to binary
    expr_matrix['Responder'] = expr_matrix[response_col].apply(
        lambda x: 1 if str(x).upper() in ['CR', 'PR', 'RESPONDER', 'RESPONSE', '1', 'YES', 'TRUE'] else 0
    )
    
    responders = expr_matrix[expr_matrix['Responder'] == 1]
    non_responders = expr_matrix[expr_matrix['Responder'] == 0]
    
    print(f"\n  Total samples with response: {len(expr_matrix)}")
    print(f"  Responders (CR/PR): {len(responders)}")
    print(f"  Non-responders (SD/PD): {len(non_responders)}")
    
    high_group = expr_matrix[expr_matrix['IIR_group'] == 'High']
    low_group = expr_matrix[expr_matrix['IIR_group'] == 'Low']
    
    high_resp = len(high_group[high_group['Responder'] == 1])
    high_non = len(high_group[high_group['Responder'] == 0])
    low_resp = len(low_group[low_group['Responder'] == 1])
    low_non = len(low_group[low_group['Responder'] == 0])
    
    print(f"\n  Contingency table:")
    print(f"                Responder  Non-responder")
    print(f"  IIR High      {high_resp:3d}        {high_non:3d}")
    print(f"  IIR Low       {low_resp:3d}        {low_non:3d}")
    
    oddsratio, p_value = fisher_exact([[high_resp, high_non], [low_resp, low_non]])
    print(f"\n  Fisher's exact test:")
    print(f"    Odds ratio: {oddsratio:.4f}")
    print(f"    p-value: {p_value:.4f}")
    
    auc = roc_auc_score(expr_matrix['Responder'], expr_matrix['IIR_score_norm'])
    print(f"\n  ROC AUC: {auc:.4f}")
    
    response_results = {
        'n_total': len(expr_matrix),
        'n_responders': len(responders),
        'n_non_responders': len(non_responders),
        'high_resp': high_resp,
        'high_non': high_non,
        'low_resp': low_resp,
        'low_non': low_non,
        'odds_ratio': oddsratio,
        'p_value': p_value,
        'auc': auc
    }

# =============================================================================
# 7. VALIDATION: SURVIVAL (PFS/OS)
# =============================================================================

print("\n" + "=" * 60)
print(" VALIDATION: SURVIVAL ANALYSIS")
print("=" * 60)

survival_results = {}

# Try to find PFS/OS columns
time_col = None
event_col = None

for col in expr_matrix.columns:
    col_lower = col.lower()
    if 'pfs' in col_lower or 'progression' in col_lower:
        if 'time' in col_lower or 'days' in col_lower or 'months' in col_lower:
            time_col = col
        elif 'event' in col_lower or 'status' in col_lower:
            event_col = col
    elif 'os' in col_lower or 'survival' in col_lower:
        if 'time' in col_lower or 'days' in col_lower or 'months' in col_lower:
            if time_col is None:
                time_col = col
        elif 'event' in col_lower or 'status' in col_lower:
            if event_col is None:
                event_col = col

print(f"  Time column: {time_col}")
print(f"  Event column: {event_col}")

if time_col and event_col and time_col in expr_matrix.columns and event_col in expr_matrix.columns:
    print("\n  Generating Kaplan-Meier curves...")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    kmf = KaplanMeierFitter()
    
    for group in ['High', 'Low']:
        mask = expr_matrix['IIR_group'] == group
        if mask.sum() > 0:
            kmf.fit(
                durations=expr_matrix.loc[mask, time_col],
                event_observed=expr_matrix.loc[mask, event_col],
                label=f'IIR {group}'
            )
            kmf.plot_survival_function(ax=ax, ci_show=True)
    
    ax.set_title('GSE91061: Survival by IIR Group', fontsize=14)
    ax.set_xlabel('Time')
    ax.set_ylabel('Survival Probability')
    ax.legend(title='IIR Group')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    fig_path = "results/figures/gse91061_km.png"
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"  ✓ KM curve saved: {fig_path}")
    
    # Log-rank test
    high_mask = expr_matrix['IIR_group'] == 'High'
    low_mask = expr_matrix['IIR_group'] == 'Low'
    
    if high_mask.sum() > 0 and low_mask.sum() > 0:
        logrank = logrank_test(
            expr_matrix.loc[high_mask, time_col],
            expr_matrix.loc[low_mask, time_col],
            expr_matrix.loc[high_mask, event_col],
            expr_matrix.loc[low_mask, event_col]
        )
        print(f"\n  Log-rank p-value: {logrank.p_value:.4f}")
        survival_results['logrank_p'] = logrank.p_value
        
        # Cox PH
        cph = CoxPHFitter()
        cph_df = expr_matrix[[time_col, event_col, 'IIR_score_norm']].dropna()
        cph.fit(cph_df, time_col, event_col)
        
        hr = np.exp(cph.params_['IIR_score_norm'])
        hr_lower = np.exp(cph.confidence_intervals_.loc['IIR_score_norm', '95% lower-bound'])
        hr_upper = np.exp(cph.confidence_intervals_.loc['IIR_score_norm', '95% upper-bound'])
        p_val = cph.summary.loc['IIR_score_norm', 'p']
        
        print(f"\n  Cox PH (IIR_score continuous):")
        print(f"    HR = {hr:.4f} (95% CI: {hr_lower:.4f} - {hr_upper:.4f})")
        print(f"    p = {p_val:.4f}")
        
        survival_results['cox_hr'] = hr
        survival_results['cox_hr_lower'] = hr_lower
        survival_results['cox_hr_upper'] = hr_upper
        survival_results['cox_p'] = p_val

# =============================================================================
# 8. SAVE RESULTS
# =============================================================================

print("\nSaving results...")

# Save IIR scores
output_file = "results/gse91061/gse91061_iir_scores.csv"
expr_matrix[['IIR_score_norm', 'IIR_group']].to_csv(output_file)
print(f"  ✓ IIR scores saved: {output_file}")

# Save summary
summary = f"""GSE91061 VALIDATION SUMMARY
========================================
Dataset: GSE91061 (Melanoma, PD-1/CTLA-4 treated)
Samples: {len(expr_matrix)}
IIR median: {median_iir:.4f}
High group: {(expr_matrix['IIR_group'] == 'High').sum()}
Low group: {(expr_matrix['IIR_group'] == 'Low').sum()}

ORR Validation:
"""
if response_results:
    summary += f"""  Odds ratio: {response_results['odds_ratio']:.4f}
  p-value: {response_results['p_value']:.4f}
  ROC AUC: {response_results['auc']:.4f}
"""
else:
    summary += "  No response data available.\n"

if survival_results:
    summary += f"""
Survival Validation:
  Log-rank p-value: {survival_results.get('logrank_p', 'N/A'):.4f}
  Cox HR: {survival_results.get('cox_hr', 'N/A'):.4f}
  Cox p-value: {survival_results.get('cox_p', 'N/A'):.4f}
"""

with open("results/gse91061/gse91061_summary.txt", "w") as f:
    f.write(summary)

print(f"  ✓ Summary saved: results/gse91061/gse91061_summary.txt")

# =============================================================================
# 9. SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print(" SUMMARY — GSE91061 VALIDATION")
print("=" * 60)

print(f"\n  Dataset: GSE91061 (Melanoma, PD-1/CTLA-4 treated)")
print(f"  Samples: {len(expr_matrix)}")

if response_results:
    print(f"\n  ORR Validation:")
    print(f"    Odds ratio: {response_results['odds_ratio']:.4f}")
    print(f"    p-value: {response_results['p_value']:.4f}")
    print(f"    ROC AUC: {response_results['auc']:.4f}")

if survival_results:
    print(f"\n  Survival Validation:")
    print(f"    Log-rank p-value: {survival_results.get('logrank_p', 'N/A'):.4f}")
    print(f"    Cox HR: {survival_results.get('cox_hr', 'N/A'):.4f}")
    print(f"    Cox p-value: {survival_results.get('cox_p', 'N/A'):.4f}")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ GSE91061 VALIDATION COMPLETE")