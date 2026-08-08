#!/usr/bin/env python3
"""
Script: 01_recalculate_tmb.py
Purpose: Recalculate TMB using standardized methods
Author: Bhaskararao Ch (Baashi)
GAP: 05 — TMB Calculation Standardization
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# 1. SETUP
# =============================================================================

os.chdir("D:/Baashi/TNBC_project/Immune_Spatial_Immunotherapy/GAP_FIXING/GAP_05_TMB_Standardization")

print("=" * 60)
print(" TMB CALCULATION STANDARDIZATION")
print(" TCGA MAF — Non-silent Variants — 38 Mb Exome")
print("=" * 60)
print(f"Started at: {pd.Timestamp.now()}\n")

os.makedirs("results", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# =============================================================================
# 2. LOAD MAF DATA
# =============================================================================

print("Loading MAF data...")

maf_file = "D:/Baashi/TNBC_project/data_raw/tcga/mutations/tcga_mutations_TNBCproxy.tsv"

if not os.path.exists(maf_file):
    print(f"ERROR: MAF file not found: {maf_file}")
    exit(1)

maf_df = pd.read_csv(maf_file, sep='\t', comment='#', low_memory=False)
print(f"  Loaded: {len(maf_df)} total variants")

# =============================================================================
# 3. FILTER TO NON-SILENT VARIANTS
# =============================================================================

print("\n" + "=" * 60)
print(" FILTERING TO NON-SILENT VARIANTS")
print("=" * 60)

# Define silent variant classifications
silent_classes = [
    'Silent', 'Intron', '3\'UTR', '5\'UTR', '3\'Flank', '5\'Flank',
    'IGR', 'RNA', 'Splice_Region'
]

non_silent = maf_df[~maf_df['Variant_Classification'].isin(silent_classes)]

print(f"  Silent variants filtered out: {len(maf_df) - len(non_silent)}")
print(f"  Non-silent variants: {len(non_silent)}")

# Non-silent breakdown
print("\n  Non-silent variant types:")
non_silent_counts = non_silent['Variant_Classification'].value_counts()
for vc, count in non_silent_counts.items():
    print(f"    {vc}: {count}")

# =============================================================================
# 4. CALCULATE TMB PER SAMPLE
# =============================================================================

print("\n" + "=" * 60)
print(" CALCULATING TMB PER SAMPLE")
print("=" * 60)

# Use Tumor_Sample_Barcode
sample_col = 'Tumor_Sample_Barcode'

if sample_col not in non_silent.columns:
    print(f"  ERROR: {sample_col} not found")
    print(f"  Available columns: {non_silent.columns[:10].tolist()}")
    exit(1)

# Count variants per sample
variant_counts = non_silent[sample_col].value_counts()
print(f"  Samples with non-silent variants: {len(variant_counts)}")

# Exome size (standard: 38 Mb for FoundationOne CDx)
EXOME_SIZE_MB = 38
print(f"  Exome size: {EXOME_SIZE_MB} Mb")

# Calculate TMB
tmb_df = pd.DataFrame({
    'sample_id': variant_counts.index,
    'n_variants': variant_counts.values,
    'TMB_mut_per_Mb': variant_counts.values / EXOME_SIZE_MB
})

print(f"\n  TMB range: {tmb_df['TMB_mut_per_Mb'].min():.2f} - {tmb_df['TMB_mut_per_Mb'].max():.2f} mut/Mb")
print(f"  TMB median: {tmb_df['TMB_mut_per_Mb'].median():.2f} mut/Mb")

# =============================================================================
# 5. EXCLUDE HYPERMUTATED SAMPLES
# =============================================================================

print("\n" + "=" * 60)
print(" EXCLUDING HYPERMUTATED SAMPLES")
print("=" * 60)

# Hypermutation threshold (common: > 20 mut/Mb)
HYPERMUT_THRESHOLD = 20

hypermut_samples = tmb_df[tmb_df['TMB_mut_per_Mb'] > HYPERMUT_THRESHOLD]
print(f"  Hypermutated samples (> {HYPERMUT_THRESHOLD} mut/Mb): {len(hypermut_samples)}")

if len(hypermut_samples) > 0:
    print(f"  Hypermutated samples and their TMB:")
    for _, row in hypermut_samples.iterrows():
        print(f"    {row['sample_id']}: {row['TMB_mut_per_Mb']:.2f} mut/Mb")
    
    # Remove hypermutated samples
    tmb_clean = tmb_df[tmb_df['TMB_mut_per_Mb'] <= HYPERMUT_THRESHOLD]
    print(f"\n  Samples after removing hypermutated: {len(tmb_clean)}")
else:
    tmb_clean = tmb_df
    print("  No hypermutated samples found.")

# =============================================================================
# 6. TMB DISTRIBUTION
# =============================================================================

print("\n" + "=" * 60)
print(" TMB DISTRIBUTION")
print("=" * 60)

print(f"\n  TMB summary (after filtering):")
print(f"    Min: {tmb_clean['TMB_mut_per_Mb'].min():.2f} mut/Mb")
print(f"    Q1: {tmb_clean['TMB_mut_per_Mb'].quantile(0.25):.2f} mut/Mb")
print(f"    Median: {tmb_clean['TMB_mut_per_Mb'].median():.2f} mut/Mb")
print(f"    Q3: {tmb_clean['TMB_mut_per_Mb'].quantile(0.75):.2f} mut/Mb")
print(f"    Max: {tmb_clean['TMB_mut_per_Mb'].max():.2f} mut/Mb")
print(f"    Samples: {len(tmb_clean)}")

# Compare to published TNBC range (median ~1.5-2.0)
print(f"\n  Published TNBC TMB range: median ~1.5-2.0 mut/Mb")
print(f"  Our TMB median: {tmb_clean['TMB_mut_per_Mb'].median():.2f} mut/Mb")

if 1.0 <= tmb_clean['TMB_mut_per_Mb'].median() <= 2.5:
    print("  ✅ TMB distribution matches published TNBC ranges!")
else:
    print("  ⚠️ TMB median outside published range. Check MAF source.")

# =============================================================================
# 7. VISUALIZE TMB DISTRIBUTION
# =============================================================================

print("\n" + "=" * 60)
print(" VISUALIZING TMB DISTRIBUTION")
print("=" * 60)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Histogram
ax = axes[0]
ax.hist(tmb_clean['TMB_mut_per_Mb'], bins=30, edgecolor='black', alpha=0.7)
ax.axvline(tmb_clean['TMB_mut_per_Mb'].median(), color='red', linestyle='--', 
           label=f'Median: {tmb_clean["TMB_mut_per_Mb"].median():.2f}')
ax.axvline(20, color='orange', linestyle='--', label='Hypermutation threshold')
ax.set_xlabel('TMB (mut/Mb)')
ax.set_ylabel('Number of Samples')
ax.set_title('TMB Distribution (TNBC, 38 Mb exome)')
ax.legend()

# Boxplot
ax = axes[1]
bp = ax.boxplot(tmb_clean['TMB_mut_per_Mb'], vert=True, patch_artist=True)
bp['boxes'][0].set_facecolor('lightblue')
ax.set_ylabel('TMB (mut/Mb)')
ax.set_title('TMB Boxplot')
ax.set_xticklabels(['TNBC'])

plt.tight_layout()
fig_path = "results/figures/tmb_distribution.png"
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
print(f"  ✓ TMB distribution plot saved: {fig_path}")

# =============================================================================
# 8. SAVE RESULTS
# =============================================================================

print("\nSaving results...")

# Save TMB data
output_file = "results/tmb_standardized.tsv"
tmb_clean.to_csv(output_file, sep='\t', index=False)
print(f"  ✓ TMB data saved: {output_file}")

# Also save original TMB for comparison
tmb_df.to_csv("results/tmb_original.tsv", sep='\t', index=False)

# Create summary
summary = f"""TMB STANDARDIZATION SUMMARY
========================================
MAF file: {maf_file}
Total variants: {len(maf_df)}
Non-silent variants: {len(non_silent)}
Exome size: {EXOME_SIZE_MB} Mb
Hypermutation threshold: {HYPERMUT_THRESHOLD} mut/Mb
Hypermutated samples removed: {len(hypermut_samples)}

TMB Distribution (n={len(tmb_clean)}):
  Min: {tmb_clean['TMB_mut_per_Mb'].min():.2f} mut/Mb
  Q1: {tmb_clean['TMB_mut_per_Mb'].quantile(0.25):.2f} mut/Mb
  Median: {tmb_clean['TMB_mut_per_Mb'].median():.2f} mut/Mb
  Q3: {tmb_clean['TMB_mut_per_Mb'].quantile(0.75):.2f} mut/Mb
  Max: {tmb_clean['TMB_mut_per_Mb'].max():.2f} mut/Mb

Published TNBC TMB range: median ~1.5-2.0 mut/Mb
Our TMB median: {tmb_clean['TMB_mut_per_Mb'].median():.2f} mut/Mb
"""

with open("results/tmb_summary.txt", "w") as f:
    f.write(summary)

print(f"  ✓ Summary saved: results/tmb_summary.txt")

# =============================================================================
# 9. SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print(" SUMMARY — TMB STANDARDIZATION")
print("=" * 60)

print(f"\n  Total variants: {len(maf_df)}")
print(f"  Non-silent variants: {len(non_silent)}")
print(f"  Exome size: {EXOME_SIZE_MB} Mb")
print(f"  Hypermutated samples removed: {len(hypermut_samples)}")
print(f"  Final samples: {len(tmb_clean)}")

print(f"\n  TMB median: {tmb_clean['TMB_mut_per_Mb'].median():.2f} mut/Mb")
print(f"  Published TNBC median: ~1.5-2.0 mut/Mb")

if 1.0 <= tmb_clean['TMB_mut_per_Mb'].median() <= 2.5:
    print("\n  ✅ TMB distribution matches published TNBC ranges!")
else:
    print("\n  ⚠️ TMB median outside published range.")

print(f"\nCompleted at: {pd.Timestamp.now()}")
print("\n✅ TMB STANDARDIZATION COMPLETE")