#!/usr/bin/env python3
"""
Script: 01_validate_iir_gse91061.py
Purpose: Validate IIR score on GSE91061 melanoma cohort (PD-1/CTLA-4 treated)
Author: Bhaskararao Ch (Baashi)
GAP: 02 — ICB-Treated Cohort Validation
"""

import os
import pandas as pd
import numpy as np
import gzip
import urllib.request
from scipy.stats import zscore, fisher_exact
from sklearn.metrics import roc_auc_score
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

os.makedirs("results", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# =============================================================================
# 2. LOAD OR DOWNLOAD GENE ANNOTATION
# =============================================================================

print("Loading gene annotation...")

annot_file = "D:/Baashi/TNBC_project/ICB_Data/GSE91061/hg19_knownGene_annot.tsv"

if not os.path.exists(annot_file):
    print("  Annotation file not found. Downloading from UCSC...")
    try:
        url = "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/knownGene.txt.gz"
        urllib.request.urlretrieve(url, annot_file + ".gz")
        # Unzip and read
        with gzip.open(annot_file + ".gz", 'rt') as f:
            annot_df = pd.read_csv(f, header=None, sep='\t')
        # Rename columns
        annot_df.columns = ['gene_id', 'chrom', 'strand', 'txStart', 'txEnd', 
                            'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 
                            'exonEnds', 'proteinID', 'alignID', 'geneSymbol']
        # Save as TSV
        annot_df.to_csv(annot_file, sep='\t', index=False)
        print(f"  ✓ Annotation saved: {annot_file}")
    except Exception as e:
        print(f"  ERROR downloading annotation: {e}")
        print("  Trying alternative method...")
        # Try a different URL
        url = "https://raw.githubusercontent.com/NCBI-Hackathons/Ensembl_Transcript_Usage/master/data/knownGene.txt"
        annot_df = pd.read_csv(url, sep='\t', header=None)
        annot_df.columns = ['gene_id', 'chrom', 'strand', 'txStart', 'txEnd', 
                            'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 
                            'exonEnds', 'proteinID', 'alignID', 'geneSymbol']
        annot_df.to_csv(annot_file, sep='\t', index=False)
        print(f"  ✓ Annotation saved: {annot_file}")
else:
    annot_df = pd.read_csv(annot_file, sep='\t')
    print(f"  Loaded annotation: {annot_df.shape[0]} genes")

# =============================================================================
# 3. LOAD EXPRESSION DATA (FPKM)
# =============================================================================

print("\nLoading GSE91061 expression data (FPKM)...")

fpkm_file = "D:/Baashi/TNBC_project/ICB_Data/GSE91061/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz"

if not os.path.exists(fpkm_file):
    print(f"ERROR: FPKM file not found: {fpkm_file}")
    exit(1)

# Load FPKM data
with gzip.open(fpkm_file, 'rt') as f:
    expr_df = pd.read_csv(f, index_col=0)

print(f"  Loaded: {expr_df.shape[0]} genes × {expr_df.shape[1]} samples")

# =============================================================================
# 4. MAP GENE IDs TO SYMBOLS
# =============================================================================

print("\nMapping gene IDs to symbols...")

# Get gene IDs from expression data (column names or index?)
# Usually the first column is gene IDs, rest are samples
# In this file, the index is gene IDs (they appear as numbers like 1, 10, 100)

gene_ids = expr_df.index.astype(str).tolist()

# Create mapping from gene_id to gene_symbol
gene_id_to_symbol = {}
for _, row in annot_df.iterrows():
    gene_id = str(row['gene_id'])
    gene_symbol = row['geneSymbol']
    gene_id_to_symbol[gene_id] = gene_symbol

# Map gene IDs to symbols
mapped_genes = []
for gene_id in gene_ids:
    if gene_id in gene_id_to_symbol:
        mapped_genes.append(gene_id_to_symbol[gene_id])
    else:
        mapped_genes.append(gene_id)  # Keep original if not found

# Assign mapped gene names as index
expr_df.index = mapped_genes

# Remove duplicates (keep first occurrence)
expr_df = expr_df[~expr_df.index.duplicated(keep='first')]

print(f"  Mapped {len(gene_ids)} genes to symbols")
print(f"  Unique genes: {expr_df.shape[0]}")

# =============================================================================
# 5. LOAD CLINICAL DATA
# =============================================================================

print("\nLoading clinical data...")

merged_file = "D:/Baashi/TNBC_project/ICB_Data/GSE91061/GSE91061_merged_data.csv"
clinical_df = pd.read_csv(merged_file, low_memory=False)

print(f"  Loaded: {clinical_df.shape[0]} rows × {clinical_df.shape[1]} columns")

# Extract clinical columns
sample_ids = clinical_df['sample_id'].tolist()

# Expression samples
expr_samples = expr_df.columns.tolist()

print(f"  Expression samples: {len(expr_samples)}")
print(f"  Clinical samples: {len(sample_ids)}")

# Find overlapping samples
overlap = set(expr_samples) & set(sample_ids)
print(f"  Overlapping samples: {len(overlap)}")

# =============================================================================
# 6. FILTER TO OVERLAPPING SAMPLES
# =============================================================================

samples_to_keep = list(overlap)
expr_df = expr_df[samples_to_keep]
expr_matrix = expr_df.T  # samples × genes

print(f"\n  Final expression matrix: {expr_matrix.shape[0]} samples × {expr_matrix.shape[1]} genes")

# =============================================================================
# 7. FIND PD1 GENES
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

gene_cols = expr_matrix.columns.tolist()

available_genes = []
for gene in pd1_genes:
    if gene in gene_cols:
        available_genes.append(gene)
    elif gene.upper() in [g.upper() for g in gene_cols]:
        # Case-insensitive match
        for col in gene_cols:
            if col.upper() == gene.upper():
                available_genes.append(col)
                break

available_genes = list(set(available_genes))

print(f"  Available genes: {len(available_genes)}/{len(pd1_genes)}")

if len(available_genes) < 5:
    print(f"  WARNING: Only {len(available_genes)} genes found.")
    print(f"  Available genes: {available_genes[:10]}")

# =============================================================================
# 8. COMPUTE IIR SCORE
# =============================================================================

print("\nComputing IIR score...")

if len(available_genes) >= 3:
    # Extract expression for available genes
    gene_expr = expr_matrix[available_genes].copy()
    
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
    
    # IIR score
    expr_matrix['IIR_score'] = (expr_matrix['PD1_signature'] + expr_matrix['ImmuneScore_norm']) / 2
    expr_matrix['IIR_score_norm'] = (expr_matrix['IIR_score'] - expr_matrix['IIR_score'].min()) / \
                                    (expr_matrix['IIR_score'].max() - expr_matrix['IIR_score'].min() + 1e-10)
    
    median_iir = expr_matrix['IIR_score_norm'].median()
    expr_matrix['IIR_group'] = np.where(expr_matrix['IIR_score_norm'] > median_iir, 'High', 'Low')
    
    print(f"  IIR range: {expr_matrix['IIR_score_norm'].min():.4f} - {expr_matrix['IIR_score_norm'].max():.4f}")
    print(f"  IIR median: {median_iir:.4f}")
    print(f"  High group: {(expr_matrix['IIR_group'] == 'High').sum()}")
    print(f"  Low group: {(expr_matrix['IIR_group'] == 'Low').sum()}")
else:
    print("  ERROR: Not enough genes found to compute IIR score.")
    print(f"  Found {len(available_genes)} genes, need at least 3.")
    exit(1)

# =============================================================================
# 9. MERGE WITH CLINICAL DATA
# =============================================================================

print("\nMerging with clinical data...")

# Create clinical mapping
clinical_dict = {}
for idx, row in clinical_df.iterrows():
    sample_id = row['sample_id']
    clinical_dict[sample_id] = {
        'response': row.get('response', None),
        'responder': row.get('responder', None),
        'visit_type': row.get('visit_type', None),
        'raw_response': row.get('raw_response', None)
    }

# Add clinical data to expression dataframe
expr_matrix['response'] = expr_matrix.index.map(
    lambda x: clinical_dict.get(x, {}).get('response', None)
)
expr_matrix['responder'] = expr_matrix.index.map(
    lambda x: clinical_dict.get(x, {}).get('responder', None)
)
expr_matrix['visit_type'] = expr_matrix.index.map(
    lambda x: clinical_dict.get(x, {}).get('visit_type', None)
)

# =============================================================================
# 10. VALIDATION: RESPONSE (ORR)
# =============================================================================

print("\n" + "=" * 60)
print(" VALIDATION: OBJECTIVE RESPONSE RATE (ORR)")
print("=" * 60)

response_results = None

if 'responder' in expr_matrix.columns:
    response_data = expr_matrix[['responder', 'IIR_group', 'IIR_score_norm', 'visit_type']].dropna()
    
    response_data['Responder'] = response_data['responder'].apply(
        lambda x: 1 if str(x).upper() in ['1', 'TRUE', 'YES', 'RESPONDER'] else 0
    )
    
    responders = response_data[response_data['Responder'] == 1]
    non_responders = response_data[response_data['Responder'] == 0]
    
    print(f"\n  Total samples with response: {len(response_data)}")
    print(f"  Responders: {len(responders)}")
    print(f"  Non-responders: {len(non_responders)}")
    
    high_group = response_data[response_data['IIR_group'] == 'High']
    low_group = response_data[response_data['IIR_group'] == 'Low']
    
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
    
    auc = roc_auc_score(response_data['Responder'], response_data['IIR_score_norm'])
    print(f"\n  ROC AUC: {auc:.4f}")
    
    response_results = {
        'dataset': 'GSE91061',
        'n_total': len(response_data),
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
# 11. SAVE RESULTS
# =============================================================================

print("\n" + "=" * 60)
print(" SAVING RESULTS")
print("=" * 60)

iir_output = "results/gse91061_iir_scores.csv"
expr_matrix[['IIR_score_norm', 'IIR_group', 'visit_type', 'response', 'responder']].to_csv(iir_output)
print(f"  ✓ IIR scores saved: {iir_output}")

if response_results:
    resp_output = "results/gse91061_response_results.csv"
    pd.DataFrame([response_results]).to_csv(resp_output, index=False)
    print(f"  ✓ Response results saved: {resp_output}")

# =============================================================================
# 12. SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print(" SUMMARY — GSE91061 VALIDATION")
print("=" * 60)

print(f"\n  Dataset: GSE91061 (Melanoma, PD-1/CTLA-4 treated)")
print(f"  Samples: {len(expr_matrix)}")
print(f"  Genes used: {len(available_genes)}")

if response_results:
    print(f"\n  ORR Validation:")
    print(f"    Odds ratio: {response_results['odds_ratio']:.4f}")
    print(f"    p-value: {response_results['p_value']:.4f}")
    print(f"    ROC AUC: {response_results['auc']:.4f}")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ GSE91061 VALIDATION COMPLETE")