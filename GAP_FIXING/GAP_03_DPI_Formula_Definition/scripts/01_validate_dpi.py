#!/usr/bin/env python3
"""
Script: 01_validate_dpi.py
Purpose: Define and validate DPI formula with biological markers
Author: Bhaskararao Ch (Baashi)
GAP: 03 — DPI Formula Definition
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# 1. SETUP
# =============================================================================

os.chdir("D:/Baashi/TNBC_project/Immune_Spatial_Immunotherapy/GAP_FIXING/GAP_03_DPI_Formula_Definition")

print("=" * 60)
print(" DPI FORMULA VALIDATION")
print(" Drug Penetration Index — Biological Validation")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

os.makedirs("results", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# =============================================================================
# 2. LOAD DATA
# =============================================================================

print("Loading data...")

# Load Spatial Metrics
spatial_file = "../../Phase_IV_Spatial_Metrics/results/spatial_metrics/spatial_metrics_with_clin_TNBC_like_TCGA.tsv"
if not os.path.exists(spatial_file):
    spatial_file = "../../../Immune_Spatial_Immunotherapy/Phase_IV_Spatial_Metrics/results/spatial_metrics/spatial_metrics_with_clin_TNBC_like_TCGA.tsv"

if not os.path.exists(spatial_file):
    print(f"ERROR: Spatial file not found: {spatial_file}")
    exit(1)

spatial_df = pd.read_csv(spatial_file, sep='\t')
print(f"  Spatial metrics: {spatial_df.shape[0]} samples")

# Load Immune Scores
immune_file = "../../Phase_I_Immune_Deconv/results/immune_scores/immune_hot_cold_TNBC_like_TCGA_OS.tsv"
if not os.path.exists(immune_file):
    immune_file = "../../../Immune_Spatial_Immunotherapy/Phase_I_Immune_Deconv/results/immune_scores/immune_hot_cold_TNBC_like_TCGA_OS.tsv"

if not os.path.exists(immune_file):
    print("ERROR: Immune scores not found.")
    exit(1)

immune_df = pd.read_csv(immune_file, sep='\t')
print(f"  Immune scores: {immune_df.shape[0]} samples")

# =============================================================================
# 3. MERGE DATA
# =============================================================================

print("\nMerging data...")

merged_df = pd.merge(spatial_df, immune_df, on='submitter_id', how='inner', suffixes=('_spatial', '_immune'))
print(f"  Merged: {merged_df.shape[0]} samples")

# =============================================================================
# 4. COMPUTE DPI
# =============================================================================

print("\n" + "=" * 60)
print(" DPI COMPUTATION")
print("=" * 60)

# Find columns
immune_col = 'ImmuneScore_immune'
stroma_col = 'Stroma_composite'
exclusion_col = 'Immune_Exclusion_index'

print(f"  Using ImmuneScore column: {immune_col}")
print(f"  Stroma column: {stroma_col}")
print(f"  Exclusion column: {exclusion_col}")

# Normalize ImmuneScore
immune_min = merged_df[immune_col].min()
immune_max = merged_df[immune_col].max()
merged_df['ImmuneScore_norm'] = (merged_df[immune_col] - immune_min) / (immune_max - immune_min + 1e-10)
print(f"  ImmuneScore_norm: {merged_df['ImmuneScore_norm'].min():.4f} - {merged_df['ImmuneScore_norm'].max():.4f}")

# Normalize Stroma
if stroma_col:
    stroma_min = merged_df[stroma_col].min()
    stroma_max = merged_df[stroma_col].max()
    merged_df['Stroma_composite_norm'] = (merged_df[stroma_col] - stroma_min) / (stroma_max - stroma_min + 1e-10)
else:
    merged_df['Stroma_composite_norm'] = 0.5
print(f"  Stroma_composite_norm: {merged_df['Stroma_composite_norm'].min():.4f} - {merged_df['Stroma_composite_norm'].max():.4f}")

# ---- DPI ORIGINAL ----
merged_df['DPI_original'] = 0.6 * (1 - merged_df[exclusion_col]) + 0.4 * (1 - merged_df['Stroma_composite_norm'])
print(f"  Original DPI: {merged_df['DPI_original'].min():.4f} - {merged_df['DPI_original'].max():.4f}")

# ---- DPI UPDATED ----
merged_df['DPI_updated'] = (merged_df['ImmuneScore_norm'] * (1 - merged_df[exclusion_col])) / \
                           (1 + merged_df['Stroma_composite_norm'])
dpi_min = merged_df['DPI_updated'].min()
dpi_max = merged_df['DPI_updated'].max()
merged_df['DPI_updated_norm'] = (merged_df['DPI_updated'] - dpi_min) / (dpi_max - dpi_min + 1e-10)
print(f"  Updated DPI: {merged_df['DPI_updated_norm'].min():.4f} - {merged_df['DPI_updated_norm'].max():.4f}")

# =============================================================================
# 5. LOAD EXPRESSION DATA
# =============================================================================

print("\n" + "=" * 60)
print(" LOADING EXPRESSION DATA")
print("=" * 60)

expr_file = "D:/Baashi/TNBC_project/ai_ml_drug_discovery/data/tcga_expr_vst.tsv"

if not os.path.exists(expr_file):
    expr_file = "D:/Baashi/TNBC_project/ai_ml_drug_discovery/data/expression/tcga_expr_vst.tsv"

if not os.path.exists(expr_file):
    print("ERROR: Expression file not found.")
    exit(1)

print(f"  Loading: {expr_file}")
expr_df = pd.read_csv(expr_file, sep='\t', index_col=0)
print(f"  Expression loaded: {expr_df.shape[0]} genes × {expr_df.shape[1]} samples")

# =============================================================================
# 6. BUILD ENSEMBL TO GENE SYMBOL MAPPING
# =============================================================================

print("\n" + "=" * 60)
print(" BUILDING GENE SYMBOL MAPPING")
print("=" * 60)

# Strip version numbers from Ensembl IDs
ensembl_ids = expr_df.index.tolist()
ensembl_stripped = [g.split('.')[0] for g in ensembl_ids]

print(f"  Stripped Ensembl IDs: {ensembl_stripped[:5]}")

# Create a mapping from Ensembl ID to Gene Symbol using a local file
# First, try to find a mapping file in your project
mapping_file = "D:/Baashi/TNBC_project/ai_ml_drug_discovery/data/gene_mapping.tsv"

if not os.path.exists(mapping_file):
    # Try alternative locations
    mapping_file = "D:/Baashi/TNBC_project/data_raw/tcga/gene_mapping.tsv"

if not os.path.exists(mapping_file):
    # Create mapping from available data
    print("  No mapping file found. Creating mapping from gene IDs...")
    
    # Build a mapping from Ensembl IDs to gene symbols
    # The gene IDs in the expression file are like "ENSG00000000003.15"
    # We'll try to find gene symbols from the gene IDs
    
    # Use the ENSEMBL IDs as is, but keep a mapping for common genes
    # We'll manually map known genes using a dictionary
    gene_mapping = {}
    
    # Common gene mappings (Ensembl ID -> Gene Symbol)
    # We'll build this from the gene IDs and known gene symbols
    known_genes = {
        'ENSG00000133639': 'VEGFA',
        'ENSG00000173511': 'VEGFB',
        'ENSG00000150630': 'VEGFC',
        'ENSG00000163462': 'VEGFD',
        'ENSG00000261371': 'PECAM1',  # CD31
        'ENSG00000100644': 'HIF1A',
        'ENSG00000110711': 'HIF1AN',
        'ENSG00000100888': 'EGLN1',
        'ENSG00000129521': 'EGLN3',
        'ENSG00000167036': 'VHL',
        'ENSG00000103126': 'SLC2A1',
        'ENSG00000162383': 'PDK1',
        'ENSG00000134333': 'LDHA',
        'ENSG00000078098': 'FAP',
        'ENSG00000136236': 'PDPN',
        'ENSG00000134853': 'PDGFRA',
        'ENSG00000113721': 'PDGFRB',
        'ENSG00000171812': 'COL1A1',
        'ENSG00000164692': 'COL1A2',
        'ENSG00000126016': 'VIM'
    }
    
    # Create mapping for all genes
    for ensembl_id in ensembl_ids:
        stripped = ensembl_id.split('.')[0]
        if stripped in known_genes:
            gene_mapping[ensembl_id] = known_genes[stripped]
        else:
            # Keep the Ensembl ID as fallback
            gene_mapping[ensembl_id] = ensembl_id
    
    print(f"  Created mapping with {len(gene_mapping)} genes")
else:
    # Load mapping file
    map_df = pd.read_csv(mapping_file, sep='\t')
    gene_mapping = dict(zip(map_df.iloc[:, 0], map_df.iloc[:, 1]))
    print(f"  Loaded mapping with {len(gene_mapping)} genes")

# Apply mapping to expression data
expr_df.index = expr_df.index.map(lambda x: gene_mapping.get(x, x))

# Show mapped gene names
print(f"  Mapped gene names (first 10): {expr_df.index[:10].tolist()}")

# =============================================================================
# 7. EXTRACT PATIENT IDS
# =============================================================================

print("\n" + "=" * 60)
print(" EXTRACTING PATIENT IDS")
print("=" * 60)

def extract_patient_id(sample_id):
    """Extract patient ID from TCGA sample ID (TCGA-XX-XXXX-XX... -> TCGA-XX-XXXX)"""
    parts = str(sample_id).split('-')
    if len(parts) >= 3:
        return '-'.join(parts[:3])
    return str(sample_id)

# Get expression patient IDs
expr_patients = {}
for sample in expr_df.columns:
    patient = extract_patient_id(sample)
    if patient not in expr_patients:
        expr_patients[patient] = []
    expr_patients[patient].append(sample)

print(f"  Expression patients: {len(expr_patients)}")

# Get merged patient IDs
merged_df['patient_id'] = merged_df['submitter_id'].apply(extract_patient_id)
print(f"  Merged patients: {merged_df['patient_id'].nunique()}")

# Find overlap
overlap = set(merged_df['patient_id']) & set(expr_patients.keys())
print(f"  Overlap: {len(overlap)} patients")

# =============================================================================
# 8. COMPUTE MARKERS
# =============================================================================

print("\n" + "=" * 60)
print(" COMPUTING MARKERS")
print("=" * 60)

# Define gene lists (gene symbols)
gene_lists = {
    'VEGF_score': ['VEGFA', 'VEGFB', 'VEGFC', 'VEGFD'],
    'CD31': ['PECAM1'],
    'Hypoxia_score': ['HIF1A', 'HIF1AN', 'EGLN1', 'EGLN3', 'VHL', 'SLC2A1', 'PDK1', 'LDHA'],
    'CAF_score': ['FAP', 'PDPN', 'PDGFRA', 'PDGFRB']
}

# Get available gene symbols
available_symbols = set(expr_df.index)

print(f"  Available gene symbols in expression: {len(available_symbols)}")

# Check which genes are available
for marker_name, genes in gene_lists.items():
    found = [g for g in genes if g in available_symbols]
    print(f"  {marker_name}: found {len(found)}/{len(genes)} genes")

# Compute markers for each patient
marker_data = {}

for patient in overlap:
    marker_data[patient] = {}
    samples = expr_patients[patient]
    
    for marker_name, genes in gene_lists.items():
        available_genes = [g for g in genes if g in available_symbols]
        if available_genes:
            # Get expression for these genes
            expr_subset = expr_df.loc[available_genes, samples]
            # Average across genes then across samples
            if not expr_subset.empty:
                marker_data[patient][marker_name] = expr_subset.mean().mean()
            else:
                marker_data[patient][marker_name] = np.nan
        else:
            marker_data[patient][marker_name] = np.nan

# Convert to DataFrame
marker_df = pd.DataFrame.from_dict(marker_data, orient='index')
marker_df.index.name = 'patient_id'
marker_df = marker_df.reset_index()

print(f"\n  Marker DataFrame: {marker_df.shape[0]} patients × {marker_df.shape[1]} markers")

# Show marker statistics
for col in marker_df.columns:
    if col != 'patient_id':
        non_na = marker_df[col].notna().sum()
        print(f"  {col}: {non_na} non-NA values")

# =============================================================================
# 9. MERGE MARKERS WITH DPI DATA
# =============================================================================

print("\n" + "=" * 60)
print(" MERGING MARKERS WITH DPI")
print("=" * 60)

merged_df = pd.merge(merged_df, marker_df, on='patient_id', how='left')
print(f"  Merged: {merged_df.shape[0]} samples")

# =============================================================================
# 10. VALIDATION RESULTS
# =============================================================================

print("\n" + "=" * 60)
print(" VALIDATION RESULTS")
print("=" * 60)

validation_results = []
marker_list = ['VEGF_score', 'CD31', 'Hypoxia_score', 'CAF_score']

for marker in marker_list:
    if marker in merged_df.columns:
        non_na = merged_df[marker].notna().sum()
        print(f"\n  {marker}: {non_na} samples with data")
        
        for dpi_name, dpi_col in [('Original', 'DPI_original'), ('Updated', 'DPI_updated_norm')]:
            data = merged_df[[dpi_col, marker]].dropna()
            if len(data) > 5:
                corr, p = pearsonr(data[dpi_col], data[marker])
                validation_results.append({
                    'dpi_version': dpi_name,
                    'marker': marker,
                    'correlation': corr,
                    'p_value': p,
                    'n': len(data)
                })
                status = "✅" if p < 0.05 else "❌"
                print(f"    {status} {dpi_name} DPI vs {marker}: r={corr:.4f}, p={p:.4f}, n={len(data)}")
            else:
                print(f"    ⚠️ {dpi_name} DPI vs {marker}: Insufficient data (n={len(data)})")

# =============================================================================
# 11. SAVE RESULTS
# =============================================================================

print("\nSaving results...")

merged_df.to_csv("results/dpi_scores_with_markers.tsv", sep='\t', index=False)
print("  ✓ DPI scores saved: results/dpi_scores_with_markers.tsv")

if validation_results:
    pd.DataFrame(validation_results).to_csv("results/dpi_validation_results.tsv", sep='\t', index=False)
    print("  ✓ Validation results saved: results/dpi_validation_results.tsv")

# =============================================================================
# 12. SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print(" SUMMARY — DPI VALIDATION")
print("=" * 60)

print("\n  DPI Formulas:")
print("    Original: DPI = 0.6*(1-Exclusion) + 0.4*(1-Stroma_norm)")
print("    Updated:  DPI = (ImmuneScore_norm × (1-Exclusion)) / (1+Stroma_norm)")

print(f"\n  DPI range (original): {merged_df['DPI_original'].min():.4f} - {merged_df['DPI_original'].max():.4f}")
print(f"  DPI range (updated):  {merged_df['DPI_updated_norm'].min():.4f} - {merged_df['DPI_updated_norm'].max():.4f}")

if validation_results:
    print("\n  Validation Results:")
    for r in validation_results:
        status = "✅ PASS" if r['p_value'] < 0.05 else "❌ FAIL"
        print(f"    {status}: {r['dpi_version']} DPI vs {r['marker']}: r={r['correlation']:.4f}, p={r['p_value']:.4f}")
else:
    print("\n  No validation results available.")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ DPI VALIDATION COMPLETE")