#!/usr/bin/env python3
"""
Script: 01_validate_spatial_metrics.py
Purpose: Validate virtual spatial metrics with known biology
Author: Bhaskararao Ch (Baashi)
GAP: 11 — Spatial/IHC Validation
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# 1. SETUP
# =============================================================================

os.chdir("D:/Baashi/TNBC_project/Immune_Spatial_Immunotherapy/GAP_FIXING/GAP_11_Spatial_IHC_Validation")

print("=" * 60)
print(" SPATIAL METRICS VALIDATION")
print(" Virtual Spatial Metrics vs Known Biology")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

os.makedirs("results", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

def extract_patient_id(sample_id):
    """Extract patient ID from TCGA sample ID (TCGA-XX-XXXX-XX... -> TCGA-XX-XXXX)"""
    parts = str(sample_id).split('-')
    if len(parts) >= 3:
        return '-'.join(parts[:3])
    return str(sample_id)

# =============================================================================
# 3. LOAD DATA
# =============================================================================

print("Loading data...")

# Load spatial metrics from Phase IV
spatial_file = "../../Phase_IV_Spatial_Metrics/results/spatial_metrics/spatial_metrics_with_clin_TNBC_like_TCGA.tsv"

if not os.path.exists(spatial_file):
    spatial_file = "../../../Immune_Spatial_Immunotherapy/Phase_IV_Spatial_Metrics/results/spatial_metrics/spatial_metrics_with_clin_TNBC_like_TCGA.tsv"

if not os.path.exists(spatial_file):
    print("ERROR: Spatial metrics file not found")
    exit(1)

spatial_df = pd.read_csv(spatial_file, sep='\t')
print(f"  Spatial metrics: {len(spatial_df)} samples")

# Load expression data
expr_file = "D:/Baashi/TNBC_project/ai_ml_drug_discovery/data/tcga_expr_vst.tsv"

if not os.path.exists(expr_file):
    expr_file = "D:/Baashi/TNBC_project/ai_ml_drug_discovery/data/expression/tcga_expr_vst.tsv"

if not os.path.exists(expr_file):
    print("ERROR: Expression file not found")
    exit(1)

expr_df = pd.read_csv(expr_file, sep='\t', index_col=0)
print(f"  Expression data: {expr_df.shape[0]} genes × {expr_df.shape[1]} samples")
print(f"  First 5 gene names: {expr_df.index[:5].tolist()}")

# =============================================================================
# 4. BUILD ENSEMBL TO GENE SYMBOL MAPPING
# =============================================================================

print("\n" + "=" * 60)
print(" BUILDING GENE SYMBOL MAPPING")
print("=" * 60)

# Manual mapping for key genes (Ensembl ID -> Gene Symbol)
gene_mapping = {
    # CD8 genes
    'ENSG00000153563': 'CD8A',
    'ENSG00000172116': 'CD8B',
    'ENSG00000145649': 'GZMA',
    'ENSG00000100453': 'GZMB',
    'ENSG00000169429': 'PRF1',
    # CD68 genes
    'ENSG00000129226': 'CD68',
    'ENSG00000177575': 'CD163',
    'ENSG00000133794': 'MRC1',
    'ENSG00000138945': 'MSR1',
    # Stromal genes
    'ENSG00000078098': 'FAP',
    'ENSG00000136236': 'PDPN',
    'ENSG00000134853': 'PDGFRA',
    'ENSG00000113721': 'PDGFRB',
    'ENSG00000108821': 'COL1A1',
    'ENSG00000164692': 'COL1A2',
    'ENSG00000126016': 'VIM',
}

# Create mapping from expression index
expr_gene_map = {}
for ensembl_id in expr_df.index:
    stripped = ensembl_id.split('.')[0]
    if stripped in gene_mapping:
        expr_gene_map[ensembl_id] = gene_mapping[stripped]

print(f"  Mapped {len(expr_gene_map)} Ensembl IDs to gene symbols")

# =============================================================================
# 5. COMPUTE KNOWN SPATIAL MARKERS
# =============================================================================

print("\n" + "=" * 60)
print(" COMPUTING KNOWN SPATIAL MARKERS")
print("=" * 60)

# Define gene lists (gene symbols)
gene_lists = {
    'CD8_score': ['CD8A', 'CD8B', 'GZMA', 'GZMB', 'PRF1'],
    'CD68_score': ['CD68', 'CD163', 'MRC1', 'MSR1'],
    'Stromal_score': ['FAP', 'PDPN', 'PDGFRA', 'PDGFRB', 'COL1A1', 'COL1A2', 'VIM']
}

# Map symbols to Ensembl IDs
ensembl_lists = {}
for marker_name, symbols in gene_lists.items():
    ensembl_lists[marker_name] = []
    for ensembl_id, symbol in expr_gene_map.items():
        if symbol in symbols:
            ensembl_lists[marker_name].append(ensembl_id)

# Compute markers
marker_data = {}

for marker_name, ensembl_ids in ensembl_lists.items():
    print(f"  {marker_name}: {len(ensembl_ids)} genes found")
    
    if len(ensembl_ids) > 0:
        # Extract expression
        expr_subset = expr_df.loc[ensembl_ids]
        # Average across genes
        marker_expr = expr_subset.mean(axis=0)
        
        # Map to dictionary
        for sample, value in marker_expr.items():
            patient = extract_patient_id(sample)
            if patient not in marker_data:
                marker_data[patient] = {}
            marker_data[patient][marker_name] = value

print(f"\n  Mapped {len(marker_data)} patients with marker data")

# =============================================================================
# 6. MERGE WITH SPATIAL DATA
# =============================================================================

print("\n" + "=" * 60)
print(" MERGING WITH SPATIAL DATA")
print("=" * 60)

# Add patient IDs to spatial data
spatial_df['patient_id'] = spatial_df['submitter_id'].apply(extract_patient_id)

# Create marker DataFrame
marker_df = pd.DataFrame.from_dict(marker_data, orient='index')
marker_df.index.name = 'patient_id'

# Merge
merged_df = spatial_df.merge(marker_df, left_on='patient_id', right_index=True, how='left')
print(f"  Merged: {len(merged_df)} samples")

# Show marker statistics
for marker in ['CD8_score', 'CD68_score', 'Stromal_score']:
    if marker in merged_df.columns:
        non_na = merged_df[marker].notna().sum()
        print(f"  {marker}: {non_na} non-NA values")

# =============================================================================
# 7. VALIDATE SPATIAL METRICS VS MARKERS
# =============================================================================

print("\n" + "=" * 60)
print(" VALIDATING SPATIAL METRICS VS MARKERS")
print("=" * 60)

validation_results = []

comparisons = [
    ('ImmuneScore', 'CD8_score', 'Spearman'),
    ('Immune_Stroma_ratio', 'CD8_score', 'Spearman'),
    ('Immune_Exclusion_index', 'CD8_score', 'Spearman'),
    ('ImmuneScore', 'Stromal_score', 'Spearman'),
    ('Immune_Stroma_ratio', 'Stromal_score', 'Spearman'),
]

for metric, marker, method in comparisons:
    if metric in merged_df.columns and marker in merged_df.columns:
        data = merged_df[[metric, marker]].dropna()
        if len(data) > 10:
            if method == 'Spearman':
                corr, p = spearmanr(data[metric], data[marker])
            else:
                corr, p = pearsonr(data[metric], data[marker])
            
            validation_results.append({
                'spatial_metric': metric,
                'marker': marker,
                'correlation': corr,
                'p_value': p,
                'n': len(data)
            })
            
            status = "✅" if abs(corr) > 0.3 and p < 0.05 else "❌"
            print(f"  {status} {metric} vs {marker}: r={corr:.4f}, p={p:.4f}, n={len(data)}")

# =============================================================================
# 8. VISUALIZE
# =============================================================================

print("\n" + "=" * 60)
print(" VISUALIZING VALIDATION")
print("=" * 60)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

plot_pairs = [
    (0, 0, 'ImmuneScore', 'CD8_score', 'ImmuneScore vs CD8+ signature'),
    (0, 1, 'Immune_Stroma_ratio', 'CD8_score', 'Immune-Stroma Ratio vs CD8+ signature'),
    (1, 0, 'Immune_Exclusion_index', 'CD8_score', 'Exclusion Index vs CD8+ signature'),
    (1, 1, 'ImmuneScore', 'Stromal_score', 'ImmuneScore vs Stromal signature')
]

for row, col, x_var, y_var, title in plot_pairs:
    ax = axes[row, col]
    if x_var in merged_df.columns and y_var in merged_df.columns:
        data = merged_df[[x_var, y_var]].dropna()
        if len(data) > 10:
            ax.scatter(data[x_var], data[y_var], alpha=0.5)
            corr, p = spearmanr(data[x_var], data[y_var])
            ax.set_xlabel(x_var)
            ax.set_ylabel(y_var)
            ax.set_title(f'{title}\nr = {corr:.4f}, p = {p:.4f}')
            ax.grid(True, alpha=0.3)

plt.tight_layout()
fig_path = "results/figures/spatial_validation.png"
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
print(f"  ✓ Validation plots saved: {fig_path}")

# =============================================================================
# 9. SAVE RESULTS
# =============================================================================

print("\nSaving results...")

merged_df.to_csv("results/spatial_validation_data.tsv", sep='\t', index=False)
print("  ✓ Data saved")

if validation_results:
    pd.DataFrame(validation_results).to_csv("results/spatial_validation_results.tsv", sep='\t', index=False)
    print("  ✓ Results saved")

# =============================================================================
# 10. SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print(" SUMMARY — SPATIAL METRICS VALIDATION")
print("=" * 60)

print(f"\n  Samples: {len(merged_df)}")

if validation_results:
    print("\n  Validation Results:")
    for result in validation_results:
        status = "✅ PASS" if abs(result['correlation']) > 0.3 and result['p_value'] < 0.05 else "❌ FAIL"
        print(f"    {status}: {result['spatial_metric']} vs {result['marker']}: r={result['correlation']:.4f}, p={result['p_value']:.4f}")
else:
    print("\n  No validation results available.")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ SPATIAL METRICS VALIDATION COMPLETE")