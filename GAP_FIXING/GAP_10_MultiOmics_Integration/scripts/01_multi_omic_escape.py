#!/usr/bin/env python3
"""
Script: 01_multi_omic_escape.py
Purpose: Integrate methylation and CNV data into escape score
Author: Bhaskararao Ch (Baashi)
GAP: 10 — Multi-Omics Integration
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# 1. SETUP
# =============================================================================

os.chdir("D:/Baashi/TNBC_project/Immune_Spatial_Immunotherapy/GAP_FIXING/GAP_10_MultiOmics_Integration")

print("=" * 60)
print(" MULTI-OMIC ESCAPE SCORE")
print(" Integrating Methylation + CNV + Expression")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

os.makedirs("results", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# =============================================================================
# 2. LOAD DATA
# =============================================================================

print("Loading data...")

# Load HLA escape data
escape_file = "../../Phase_X_HLA_ImmuneEscape/results/hla/HLA_escape_TNBC_like_TCGA.tsv"
if not os.path.exists(escape_file):
    escape_file = "../../../Immune_Spatial_Immunotherapy/Phase_X_HLA_ImmuneEscape/results/hla/HLA_escape_TNBC_like_TCGA.tsv"

if not os.path.exists(escape_file):
    print("ERROR: HLA escape file not found")
    exit(1)

escape_df = pd.read_csv(escape_file, sep='\t')
print(f"  HLA escape data: {len(escape_df)} samples")

# Load CNV data (GISTIC)
cnv_file = "data/TCGA_BRCA_GISTIC_thresholded.csv"
if not os.path.exists(cnv_file):
    print("ERROR: GISTIC file not found")
    exit(1)

cnv_df = pd.read_csv(cnv_file, index_col=0)
print(f"  CNV data: {cnv_df.shape[0]} genes × {cnv_df.shape[1]} samples")

# Load methylation data - FIXED
methyl_file = "D:/Baashi/TNBC_project/data_raw/tcga/methylation/processed/methylation_beta_TNBCproxy.tsv"
methyl_dict = {}

if os.path.exists(methyl_file):
    # Read with probe_id as index
    methyl_df = pd.read_csv(methyl_file, sep='\t', index_col=0)
    print(f"  Methylation data: {methyl_df.shape[0]} probes × {methyl_df.shape[1]} samples")
    print(f"  First 5 probe IDs: {methyl_df.index[:5].tolist()}")
    
    # MHC-I keywords to search for
    mhc1_keywords = ['HLA-A', 'HLA-B', 'HLA-C', 'B2M', 'TAP1', 'TAP2']
    
    # Find probes containing these keywords (case-insensitive)
    mhc1_probes = []
    for probe in methyl_df.index:
        probe_str = str(probe).upper()
        for keyword in mhc1_keywords:
            if keyword.upper() in probe_str:
                mhc1_probes.append(probe)
                break
    
    print(f"  MHC-I probes found: {len(mhc1_probes)}")
    
    if len(mhc1_probes) > 0:
        # Average methylation for MHC-I probes per sample
        mhc1_methyl = methyl_df.loc[mhc1_probes].mean(axis=0)
        methyl_dict = mhc1_methyl.to_dict()
        print(f"  MHC-I methylation range: {min(methyl_dict.values()):.4f} - {max(methyl_dict.values()):.4f}")
        print(f"  Methylation samples: {len(methyl_dict)}")
    else:
        print("  No MHC-I probes found. Trying partial matching...")
        # Try partial matching with any probe containing HLA or B2M
        partial_probes = []
        for probe in methyl_df.index:
            probe_str = str(probe).upper()
            if 'HLA' in probe_str or 'B2M' in probe_str or 'TAP' in probe_str:
                partial_probes.append(probe)
        
        print(f"  Partial match probes found: {len(partial_probes)}")
        if len(partial_probes) > 0:
            mhc1_methyl = methyl_df.loc[partial_probes].mean(axis=0)
            methyl_dict = mhc1_methyl.to_dict()
            print(f"  MHC-I methylation range: {min(methyl_dict.values()):.4f} - {max(methyl_dict.values()):.4f}")
else:
    print("  Methylation data not found")

# =============================================================================
# 3. EXTRACT IMMUNE-RELATED CNVs
# =============================================================================

print("\n" + "=" * 60)
print(" EXTRACTING IMMUNE-RELATED CNVs")
print("=" * 60)

immune_cnv_genes = [
    'IFNA1', 'IFNA2', 'IFNA4', 'IFNA5', 'IFNA6', 'IFNA7', 'IFNA8', 'IFNA10', 'IFNA13', 'IFNA14', 'IFNA16', 'IFNA17', 'IFNA21',
    'HLA-A', 'HLA-B', 'HLA-C', 'B2M', 'TAP1', 'TAP2',
    'HLA-DRA', 'HLA-DRB1', 'HLA-DQB1', 'HLA-DQA1', 'HLA-DPA1', 'HLA-DPB1',
    'PDCD1', 'CD274', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT', 'BTLA'
]

available_cnv_genes = [g for g in immune_cnv_genes if g in cnv_df.index]
print(f"  Available immune-related genes in CNV: {len(available_cnv_genes)}")

if len(available_cnv_genes) > 0:
    cnv_subset = cnv_df.loc[available_cnv_genes]
    cnv_loss = (cnv_subset < 0).sum(axis=0)
    cnv_gain = (cnv_subset > 0).sum(axis=0)
    cnv_total = (cnv_subset != 0).sum(axis=0)
    
    print(f"  CNV loss events: {cnv_loss.min()} - {cnv_loss.max()}")
    print(f"  CNV gain events: {cnv_gain.min()} - {cnv_gain.max()}")
else:
    cnv_subset = pd.DataFrame()
    cnv_loss = pd.Series()
    cnv_gain = pd.Series()
    cnv_total = pd.Series()

# =============================================================================
# 4. MAP DATA TO PATIENT IDs
# =============================================================================

print("\n" + "=" * 60)
print(" MAPPING DATA TO PATIENT IDs")
print("=" * 60)

def extract_patient_id(sample_id):
    parts = str(sample_id).split('-')
    if len(parts) >= 3:
        return '-'.join(parts[:3])
    return str(sample_id)

multi_df = escape_df.copy()
multi_df['patient_id'] = multi_df['submitter_id'].apply(extract_patient_id)

# Map CNV data
cnv_dict = {}
if len(cnv_subset) > 0:
    for sample in cnv_subset.columns:
        patient = extract_patient_id(sample)
        cnv_dict[patient] = {
            'cnv_loss': cnv_loss[sample],
            'cnv_gain': cnv_gain[sample],
            'cnv_total': cnv_total[sample]
        }
    
    multi_df['cnv_loss'] = multi_df['patient_id'].map(lambda x: cnv_dict.get(x, {}).get('cnv_loss', np.nan))
    multi_df['cnv_gain'] = multi_df['patient_id'].map(lambda x: cnv_dict.get(x, {}).get('cnv_gain', np.nan))
    multi_df['cnv_total'] = multi_df['patient_id'].map(lambda x: cnv_dict.get(x, {}).get('cnv_total', np.nan))
    print(f"  CNV data mapped: {multi_df['cnv_total'].notna().sum()} patients")

# Map methylation data
if methyl_dict:
    methyl_patients = {}
    for sample, value in methyl_dict.items():
        patient = extract_patient_id(sample)
        if patient not in methyl_patients:
            methyl_patients[patient] = []
        methyl_patients[patient].append(value)
    
    methyl_avg = {p: np.mean(vals) for p, vals in methyl_patients.items() if vals}
    multi_df['methylation_score'] = multi_df['patient_id'].map(methyl_avg)
    print(f"  Methylation data mapped: {multi_df['methylation_score'].notna().sum()} patients")

# =============================================================================
# 5. BUILD MULTI-OMIC ESCAPE SCORE
# =============================================================================

print("\n" + "=" * 60)
print(" BUILDING MULTI-OMIC ESCAPE SCORE")
print("=" * 60)

# Component 1: MHC-I expression
if 'MHC_I_axis' in multi_df.columns:
    mhc_min = multi_df['MHC_I_axis'].min()
    mhc_max = multi_df['MHC_I_axis'].max()
    multi_df['MHC_I_expr_norm'] = (multi_df['MHC_I_axis'] - mhc_min) / (mhc_max - mhc_min + 1e-10)
    multi_df['MHC_I_escape_score'] = 1 - multi_df['MHC_I_expr_norm']
    print("  Component 1: MHC-I expression (transcriptomic)")

# Component 2: Methylation
if 'methylation_score' in multi_df.columns:
    meth_min = multi_df['methylation_score'].min()
    meth_max = multi_df['methylation_score'].max()
    multi_df['methylation_escape_score'] = (multi_df['methylation_score'] - meth_min) / (meth_max - meth_min + 1e-10)
    print("  Component 2: MHC-I methylation (epigenetic)")

# Component 3: CNV loss
if 'cnv_loss' in multi_df.columns:
    cnv_min = multi_df['cnv_loss'].min()
    cnv_max = multi_df['cnv_loss'].max()
    multi_df['cnv_escape_score'] = (multi_df['cnv_loss'] - cnv_min) / (cnv_max - cnv_min + 1e-10)
    print("  Component 3: CNV loss (genomic)")

# Build multi-omic escape score
weights = {'MHC_I_escape_score': 0.4, 'methylation_escape_score': 0.3, 'cnv_escape_score': 0.3}
available_components = [c for c in weights.keys() if c in multi_df.columns]

if available_components:
    total_weight = sum(weights[c] for c in available_components)
    multi_df['multi_omic_escape_score'] = sum(
        multi_df[c] * weights[c] / total_weight for c in available_components
    )
    print(f"\n  Multi-omic escape score built from {len(available_components)} components")
    print(f"    Range: {multi_df['multi_omic_escape_score'].min():.4f} - {multi_df['multi_omic_escape_score'].max():.4f}")
else:
    print("  No components available")

# =============================================================================
# 6. COMPARE TO MHC-I ALONE
# =============================================================================

print("\n" + "=" * 60)
print(" COMPARING TO MHC-I ALONE")
print("=" * 60)

if 'multi_omic_escape_score' in multi_df.columns and 'MHC_I_escape_score' in multi_df.columns:
    mhc_median = multi_df['MHC_I_escape_score'].median()
    multi_median = multi_df['multi_omic_escape_score'].median()
    
    multi_df['mhc_i_escape_group'] = np.where(
        multi_df['MHC_I_escape_score'] > mhc_median, 'Escape', 'No_escape'
    )
    multi_df['multi_omic_escape_group'] = np.where(
        multi_df['multi_omic_escape_score'] > multi_median, 'Escape', 'No_escape'
    )
    
    print(f"\n  MHC-I alone: {(multi_df['mhc_i_escape_group'] == 'Escape').sum()} escape")
    print(f"  Multi-omic: {(multi_df['multi_omic_escape_group'] == 'Escape').sum()} escape")
    
    additional = (multi_df['multi_omic_escape_group'] == 'Escape') & (multi_df['mhc_i_escape_group'] == 'No_escape')
    print(f"\n  Additional escape patients identified by multi-omic: {additional.sum()}")

# =============================================================================
# 7. SAVE RESULTS
# =============================================================================

print("\nSaving results...")
multi_df.to_csv("results/multi_omic_escape_data.tsv", sep='\t', index=False)
print("  ✓ Data saved")

summary = f"""MULTI-OMIC ESCAPE SCORE SUMMARY
========================================
Samples: {len(multi_df)}

Escape Groups:
  MHC-I alone: {(multi_df['mhc_i_escape_group'] == 'Escape').sum() if 'mhc_i_escape_group' in multi_df.columns else 0}
  Multi-omic: {(multi_df['multi_omic_escape_group'] == 'Escape').sum() if 'multi_omic_escape_group' in multi_df.columns else 0}
"""

with open("results/multi_omic_summary.txt", "w") as f:
    f.write(summary)
print("  ✓ Summary saved")

print("\n✅ MULTI-OMIC ESCAPE SCORE COMPLETE")