#!/usr/bin/env python3
"""
Script: 01_decision_curve_analysis.py
Purpose: Decision curve analysis for IIR vs PD-L1 vs TMB
Author: Bhaskararao Ch (Baashi)
GAP: 07 — Clinical Decision Curves
"""

import os
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# 1. SETUP
# =============================================================================

os.chdir("D:/Baashi/TNBC_project/Immune_Spatial_Immunotherapy/GAP_FIXING/GAP_07_Decision_Curves")

print("=" * 60)
print(" DECISION CURVE ANALYSIS")
print(" Net Benefit — IIR vs PD-L1 vs TMB")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

os.makedirs("results", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# =============================================================================
# 2. LOAD DATA
# =============================================================================

print("Loading data...")

data_file = "../GAP_06_Biomarker_Comparison/results/biomarker_data.tsv"

if not os.path.exists(data_file):
    data_file = "../../GAP_06_Biomarker_Comparison/results/biomarker_data.tsv"

if not os.path.exists(data_file):
    print(f"ERROR: Data file not found: {data_file}")
    exit(1)

df = pd.read_csv(data_file, sep='\t')
print(f"  Loaded: {len(df)} samples")

# =============================================================================
# 3. PREPARE DATA FOR DCA
# =============================================================================

print("\nPreparing data for DCA...")

# Define binary outcome: 5-year mortality (OS_MONTHS < 60)
df['outcome_5yr'] = np.where((df['OS_MONTHS'] < 60) & (df['OS_STATUS'] == 1), 1, 0)
df['outcome_5yr'] = df['outcome_5yr'].fillna(0).astype(int)

print(f"  5-year mortality: {df['outcome_5yr'].sum()} events ({df['outcome_5yr'].mean():.1%})")

# =============================================================================
# 4. DCA IMPLEMENTATION
# =============================================================================

def calculate_net_benefit(y_true, y_pred, threshold):
    """Calculate net benefit at a given threshold."""
    y_pred_binary = (y_pred >= threshold).astype(int)
    tp = np.sum((y_true == 1) & (y_pred_binary == 1))
    fp = np.sum((y_true == 0) & (y_pred_binary == 1))
    n = len(y_true)
    net_benefit = (tp / n) - (fp / n) * (threshold / (1 - threshold))
    return net_benefit

def calculate_decision_curve(y_true, y_pred, thresholds):
    """Calculate net benefit across thresholds."""
    net_benefits = []
    for threshold in thresholds:
        nb = calculate_net_benefit(y_true, y_pred, threshold)
        net_benefits.append(nb)
    return np.array(net_benefits)

def calculate_auc_nb(thresholds, net_benefit):
    """Calculate area under the net benefit curve."""
    max_nb = max(net_benefit)
    min_nb = min(net_benefit)
    if max_nb > min_nb:
        normalized = (net_benefit - min_nb) / (max_nb - min_nb)
    else:
        normalized = net_benefit * 0
    # Use np.trapezoid (new) or np.trapz (old)
    try:
        return np.trapezoid(normalized, thresholds)
    except AttributeError:
        return np.trapz(normalized, thresholds)

# =============================================================================
# 5. TRAIN MODELS FOR DCA
# =============================================================================

print("\n" + "=" * 60)
print(" TRAINING MODELS FOR DCA")
print("=" * 60)

predictors = {
    'IIR': 'IIR_score_norm',
    'PD-L1': 'PD1_PDL1_signature',
    'ImmuneScore': 'ImmuneScore_norm',
    'TMB': 'TMB'
}

dca_df = df[['outcome_5yr'] + [v for v in predictors.values() if v in df.columns]].dropna()
print(f"\n  DCA data: {len(dca_df)} samples")

models = {}
for name, var in predictors.items():
    if var in dca_df.columns:
        X = dca_df[[var]].values.reshape(-1, 1)
        y = dca_df['outcome_5yr'].values
        model = LogisticRegression(max_iter=1000)
        model.fit(X, y)
        models[name] = model
        print(f"  {name}: trained (coef={model.coef_[0][0]:.4f})")

# =============================================================================
# 6. CALCULATE DECISION CURVES
# =============================================================================

print("\n" + "=" * 60)
print(" CALCULATING DECISION CURVES")
print("=" * 60)

thresholds = np.linspace(0.05, 0.50, 46)
results = {}

for name, model in models.items():
    X = dca_df[[predictors[name]]].values.reshape(-1, 1)
    y_pred = model.predict_proba(X)[:, 1]
    net_benefit = calculate_decision_curve(dca_df['outcome_5yr'].values, y_pred, thresholds)
    results[name] = net_benefit
    print(f"  {name}: net benefit calculated")

# Treat all
treat_all_nb = []
for threshold in thresholds:
    tp_all = np.sum(dca_df['outcome_5yr'] == 1)
    fp_all = len(dca_df)
    n = len(dca_df)
    nb = (tp_all / n) - (fp_all / n) * (threshold / (1 - threshold))
    treat_all_nb.append(nb)

# Treat none
treat_none_nb = np.zeros(len(thresholds))

# =============================================================================
# 7. VISUALIZE DECISION CURVES
# =============================================================================

print("\n" + "=" * 60)
print(" VISUALIZING DECISION CURVES")
print("=" * 60)

fig, ax = plt.subplots(figsize=(10, 6))

colors = {'IIR': '#e74c3c', 'PD-L1': '#3498db', 'ImmuneScore': '#2ecc71', 'TMB': '#f39c12'}

for name, net_benefit in results.items():
    ax.plot(thresholds * 100, net_benefit, label=name, color=colors.get(name, 'gray'), linewidth=2)

ax.plot(thresholds * 100, treat_all_nb, label='Treat All', color='black', linestyle='--', linewidth=1)
ax.plot(thresholds * 100, treat_none_nb, label='Treat None', color='black', linestyle=':', linewidth=1)

ax.set_xlabel('Threshold Probability (%)', fontsize=14)
ax.set_ylabel('Net Benefit', fontsize=14)
ax.set_title('Decision Curve Analysis: IIR vs PD-L1 vs TMB', fontsize=16)
ax.legend(loc='best')
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 50)
ax.set_ylim(-0.1, 0.4)

plt.tight_layout()
fig_path = "results/figures/decision_curves.png"
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
print(f"  ✓ Decision curves saved: {fig_path}")

# =============================================================================
# 8. NET BENEFIT COMPARISON
# =============================================================================

print("\n" + "=" * 60)
print(" NET BENEFIT COMPARISON")
print("=" * 60)

key_thresholds = [0.10, 0.20, 0.30, 0.40, 0.50]
key_indices = [np.argmin(np.abs(thresholds - t)) for t in key_thresholds]

print("\n  Net benefit at key thresholds:")
print(f"  {'Threshold':<12} {'IIR':<10} {'PD-L1':<10} {'TMB':<10} {'Treat All':<10}")
print("  " + "-" * 55)

for t, idx in zip(key_thresholds, key_indices):
    iir_nb = results.get('IIR', [0])[idx] if 'IIR' in results else 0
    pd1_nb = results.get('PD-L1', [0])[idx] if 'PD-L1' in results else 0
    tmb_nb = results.get('TMB', [0])[idx] if 'TMB' in results else 0
    all_nb = treat_all_nb[idx]
    print(f"  {t*100:>6.0f}%    {iir_nb:>8.4f}  {pd1_nb:>8.4f}  {tmb_nb:>8.4f}  {all_nb:>8.4f}")

# =============================================================================
# 9. DCA SUMMARY STATISTICS
# =============================================================================

print("\n" + "=" * 60)
print(" DCA SUMMARY STATISTICS")
print("=" * 60)

for name, net_benefit in results.items():
    auc_nb = calculate_auc_nb(thresholds, net_benefit)
    print(f"  {name}: AUC-NB = {auc_nb:.4f}")

# =============================================================================
# 10. SAVE RESULTS
# =============================================================================

print("\nSaving results...")

for name, net_benefit in results.items():
    nb_df = pd.DataFrame({
        'threshold': thresholds * 100,
        'net_benefit': net_benefit
    })
    nb_df.to_csv(f"results/net_benefit_{name}.tsv", sep='\t', index=False)

print(f"  ✓ Net benefit data saved")

# Summary
summary = f"""DECISION CURVE ANALYSIS SUMMARY
========================================
Dataset: METABRIC TNBC (n={len(dca_df)})
Outcome: 5-year mortality ({dca_df['outcome_5yr'].sum()} events)

Net benefit at key thresholds:
"""
for t, idx in zip(key_thresholds, key_indices):
    iir_nb = results.get('IIR', [0])[idx] if 'IIR' in results else 0
    pd1_nb = results.get('PD-L1', [0])[idx] if 'PD-L1' in results else 0
    tmb_nb = results.get('TMB', [0])[idx] if 'TMB' in results else 0
    all_nb = treat_all_nb[idx]
    summary += f"  {t*100:>6.0f}%: IIR={iir_nb:.4f}, PD-L1={pd1_nb:.4f}, TMB={tmb_nb:.4f}, Treat All={all_nb:.4f}\n"

summary += f"\nAUC-NB:\n"
for name, net_benefit in results.items():
    auc_nb = calculate_auc_nb(thresholds, net_benefit)
    summary += f"  {name}: {auc_nb:.4f}\n"

with open("results/dca_summary.txt", "w") as f:
    f.write(summary)

print(f"  ✓ Summary saved: results/dca_summary.txt")

# =============================================================================
# 11. SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print(" SUMMARY — DECISION CURVE ANALYSIS")
print("=" * 60)

print(f"\n  Dataset: {len(dca_df)} samples")
print(f"  Outcome: 5-year mortality ({dca_df['outcome_5yr'].sum()} events)")

if 'IIR' in results:
    iir_nb = results['IIR']
    positive_range = np.sum(iir_nb > 0) / len(iir_nb)
    print(f"\n  IIR positive net benefit: {positive_range:.1%} of thresholds")

    if 'PD-L1' in results:
        pd1_nb = results['PD-L1']
        iir_better = np.sum(iir_nb > pd1_nb) / len(iir_nb)
        print(f"  IIR > PD-L1: {iir_better:.1%} of thresholds")

    iir_vs_all = np.sum(iir_nb > treat_all_nb) / len(iir_nb)
    print(f"  IIR > Treat All: {iir_vs_all:.1%} of thresholds")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ DECISION CURVE ANALYSIS COMPLETE")