#!/usr/bin/env python3
"""
Script: 01_compute_iir_depmap.py
Purpose: Compute IIR score on DepMap cell lines for mechanistic validation
Author: Bhaskararao Ch (Baashi)
GAP: 02 — ICB-Treated Cohort Validation (Mechanistic)
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import zscore
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# 1. SETUP
# =============================================================================

os.chdir("D:/Baashi/TNBC_project/Immune_Spatial_Immunotherapy/GAP_FIXING/GAP_02_ICB_Treated_Cohort_IMpassion")

print("=" * 60)
print(" DEPMAP IIR SCORE COMPUTATION")
print(" Mechanistic Validation — Cell Lines")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

os.makedirs("results/depmap", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# =============================================================================
# 2. LOAD EXPRESSION DATA (CSV format)
# =============================================================================

print("Loading DepMap expression data (CSV)...")

expr_file = "D:/Baashi/TNBC_project/M10_cellline_validation/inputs/DepMap_expression_raw.csv"

if not os.path.exists(expr_file):
    print(f"ERROR: Expression file not found: {expr_file}")
    exit(1)

# Load CSV with first column as index (cell line IDs)
expr_df = pd.read_csv(expr_file, index_col=0)
print(f"  Loaded: {expr_df.shape[0]} cell lines × {expr_df.shape[1]} genes")

# Get cell line IDs and gene names
cell_lines = expr_df.index.tolist()
gene_cols = expr_df.columns.tolist()

print(f"  Cell lines: {len(cell_lines)}")
print(f"  Genes: {len(gene_cols)}")
print(f"  First 5 genes: {gene_cols[:5]}")

# =============================================================================
# 3. LOAD CELL LINE METADATA
# =============================================================================

print("\nLoading cell line metadata...")

meta_file = "D:/Baashi/TNBC_project/Biomarker_Verification/inputs/external_cohorts/depmap/depmap_model.csv"

if os.path.exists(meta_file):
    meta_df = pd.read_csv(meta_file)
    print(f"  Loaded: {meta_df.shape[0]} cell lines")
    
    # Check if we can map cancer types
    if 'disease_type' in meta_df.columns:
        print("\n  Cancer types in DepMap:")
        counts = meta_df['disease_type'].value_counts().head(10)
        for ct, count in counts.items():
            print(f"    {ct}: {count}")
        
        # Filter to breast cancer if needed
        # breast_cells = meta_df[meta_df['disease_type'].str.contains('breast', case=False)]['depmap_id'].tolist()
        # print(f"  Breast cancer cell lines: {len(breast_cells)}")
else:
    meta_df = None
    print("  Metadata not found — proceeding without filtering.")

# =============================================================================
# 4. EXTRACT GENE SYMBOLS FROM COLUMN NAMES
# =============================================================================

print("\nExtracting gene symbols from column names...")

# Column names are like "TSPAN6 (7105)" — extract gene symbol
def extract_gene_symbol(col_name):
    """Extract gene symbol from 'TSPAN6 (7105)' format"""
    if '(' in col_name:
        return col_name.split('(')[0].strip()
    return col_name

gene_symbols = [extract_gene_symbol(col) for col in gene_cols]

# Create a mapping from gene symbol to column name
gene_to_col = dict(zip(gene_symbols, gene_cols))

print(f"  Extracted {len(gene_symbols)} gene symbols")
print(f"  First 10: {gene_symbols[:10]}")

# =============================================================================
# 5. FIND PD1 SIGNATURE GENES
# =============================================================================

print("\nSearching for PD1 signature genes...")

pd1_genes = [
    'PDCD1', 'CD274', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT', 'BTLA',
    'CXCL9', 'CXCL10', 'CXCL11', 'CXCL13',
    'IFNG', 'STAT1', 'IRF1',
    'CD8A', 'CD8B', 'GZMA', 'GZMB', 'PRF1', 'GNLY', 'NKG7',
    'CD3D', 'CD3E', 'CD3G', 'CD4',
    'HLA-DRA', 'HLA-DRB1', 'HLA-DQB1'
]

# Find available genes
available_genes = []
available_cols = []

for gene in pd1_genes:
    if gene in gene_to_col:
        available_genes.append(gene)
        available_cols.append(gene_to_col[gene])

print(f"  Available genes: {len(available_genes)}/{len(pd1_genes)}")

if len(available_genes) < 5:
    print(f"  WARNING: Only {len(available_genes)} genes found.")
    print(f"  Available: {available_genes[:10]}")

if not available_genes:
    print("  ERROR: No PD1 signature genes found.")
    print("  Check gene naming format in DepMap.")
    print(f"  Example column names: {gene_cols[:5]}")
    exit(1)

# =============================================================================
# 6. COMPUTE IIR SCORE
# =============================================================================

print("\nComputing IIR score on cell lines...")

# Extract expression for available genes
gene_expr = expr_df[available_cols]
gene_expr.columns = available_genes  # Rename to gene symbols

print(f"  Using {gene_expr.shape[1]} genes for IIR computation")

# Z-score normalize each gene across cell lines
gene_z = gene_expr.apply(zscore, axis=0)

# PD1 signature (mean of z-scores)
expr_df['PD1_signature'] = gene_z.mean(axis=1)

# ImmuneScore (cytotoxic signature)
cytotoxic_genes = ['CD8A', 'CD8B', 'GZMA', 'GZMB', 'PRF1']
avail_cytotoxic = [g for g in cytotoxic_genes if g in available_genes]

if avail_cytotoxic:
    cyt_cols = [gene_to_col[g] for g in avail_cytotoxic]
    expr_df['ImmuneScore'] = expr_df[cyt_cols].mean(axis=1)
    expr_df['ImmuneScore_norm'] = (expr_df['ImmuneScore'] - expr_df['ImmuneScore'].min()) / \
                                  (expr_df['ImmuneScore'].max() - expr_df['ImmuneScore'].min() + 1e-10)
else:
    expr_df['ImmuneScore_norm'] = expr_df['PD1_signature']

# Spatial metrics (stromal signature)
stromal_genes = ['FAP', 'PDPN', 'PDGFRA', 'PDGFRB', 'COL1A1', 'COL1A2', 'VIM']
avail_stromal = [g for g in stromal_genes if g in gene_to_col]

if avail_stromal:
    strom_cols = [gene_to_col[g] for g in avail_stromal]
    stroma_expr = expr_df[strom_cols]
    stroma_z = stroma_expr.apply(zscore, axis=0)
    expr_df['Stroma_score'] = stroma_z.mean(axis=1)
    
    total = expr_df['ImmuneScore_norm'] + expr_df['Stroma_score']
    expr_df['Immune_Stroma_ratio'] = np.where(
        total > 0,
        expr_df['ImmuneScore_norm'] / total,
        0.5
    )
    expr_df['Immune_Exclusion_index'] = 1 - expr_df['Immune_Stroma_ratio']
else:
    expr_df['Immune_Stroma_ratio'] = 0.5
    expr_df['Immune_Exclusion_index'] = 0.5

# DPI
expr_df['DPI'] = (expr_df['ImmuneScore_norm'] * expr_df['Immune_Stroma_ratio']) / \
                 (1 + expr_df['Immune_Exclusion_index'])
expr_df['DPI_norm'] = (expr_df['DPI'] - expr_df['DPI'].min()) / \
                      (expr_df['DPI'].max() - expr_df['DPI'].min() + 1e-10)

# IIR Score
components = ['ImmuneScore_norm', 'PD1_signature', 'Immune_Exclusion_index', 'DPI_norm']

for comp in components:
    expr_df[f'{comp}_z'] = zscore(expr_df[comp])

expr_df['Exclusion_inv_z'] = -zscore(expr_df['Immune_Exclusion_index'])

expr_df['IIR_score'] = (expr_df['ImmuneScore_norm_z'] + 
                        expr_df['PD1_signature_z'] +
                        expr_df['Exclusion_inv_z'] +
                        expr_df['DPI_norm_z']) / 4

expr_df['IIR_score_norm'] = (expr_df['IIR_score'] - expr_df['IIR_score'].min()) / \
                            (expr_df['IIR_score'].max() - expr_df['IIR_score'].min() + 1e-10)

tertiles = np.percentile(expr_df['IIR_score_norm'], [33.33, 66.67])
expr_df['IIR_group'] = pd.cut(
    expr_df['IIR_score_norm'],
    bins=[-np.inf, tertiles[0], tertiles[1], np.inf],
    labels=['Low', 'Mid', 'High']
)

print(f"  IIR_score range: {expr_df['IIR_score_norm'].min():.4f} - {expr_df['IIR_score_norm'].max():.4f}")
print(f"  IIR group distribution:")
print(expr_df['IIR_group'].value_counts().sort_index())

# =============================================================================
# 7. SAVE RESULTS
# =============================================================================

print("\nSaving results...")

# Add cell line IDs to output
result_df = expr_df[['IIR_score_norm', 'IIR_group', 'ImmuneScore_norm', 'PD1_signature', 'DPI_norm']].copy()
result_df.index.name = 'cell_line_id'
result_df = result_df.reset_index()

output_file = "results/depmap/depmap_iir_scores.csv"
result_df.to_csv(output_file, index=False)
print(f"  ✓ IIR scores saved: {output_file}")

# Summary
summary = f"""DEPMAP IIR SCORE SUMMARY
========================================
Cell lines: {len(expr_df)}
Genes used: {len(available_genes)}/{len(pd1_genes)}

IIR_score distribution:
  Min: {expr_df['IIR_score_norm'].min():.4f}
  Q1: {expr_df['IIR_score_norm'].quantile(0.25):.4f}
  Median: {expr_df['IIR_score_norm'].median():.4f}
  Q3: {expr_df['IIR_score_norm'].quantile(0.75):.4f}
  Max: {expr_df['IIR_score_norm'].max():.4f}

IIR group distribution:
  Low: {(expr_df['IIR_group'] == 'Low').sum()}
  Mid: {(expr_df['IIR_group'] == 'Mid').sum()}
  High: {(expr_df['IIR_group'] == 'High').sum()}

Component correlations with IIR_score:
"""
for comp in ['ImmuneScore_norm', 'PD1_signature', 'Immune_Exclusion_index', 'DPI_norm']:
    if comp in expr_df.columns:
        corr = expr_df['IIR_score_norm'].corr(expr_df[comp])
        summary += f"  {comp}: {corr:.4f}\n"

with open("results/depmap/depmap_iir_summary.txt", "w") as f:
    f.write(summary)

print(f"  ✓ Summary saved: results/depmap/depmap_iir_summary.txt")

print("\n" + "=" * 60)
print(" SUMMARY — DEPMAP IIR COMPUTATION")
print("=" * 60)

print(f"\n  Cell lines: {len(expr_df)}")
print(f"  Genes available: {len(available_genes)}/{len(pd1_genes)}")
print(f"  IIR_score range: {expr_df['IIR_score_norm'].min():.4f} - {expr_df['IIR_score_norm'].max():.4f}")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ DEPMAP IIR SCORE COMPUTATION COMPLETE")
print("   Next: 02_validate_iir_depmap.py")