import pandas as pd
import numpy as np

print("🧹 DATA CLEANING SCRIPT")
print("="*60)

# Load data
df = pd.read_csv("results/immune_subtypes/immune_subtypes_TNBC_like_TCGA_core.tsv", sep='\t')
print(f"📥 Loaded: {len(df)} rows, {df['submitter_id'].nunique()} unique patients")

# 1. Find duplicates
dup_mask = df.duplicated(subset=['submitter_id'], keep=False)
dup_count = dup_mask.sum()
dup_patients = df[dup_mask]['submitter_id'].nunique()
print(f"\n🔍 DUPLICATES FOUND: {dup_count} duplicate rows for {dup_patients} patients")

if dup_count > 0:
    print("\n📋 Sample duplicates (showing difference in PD-L1 signature):")
    for pid in df[dup_mask]['submitter_id'].unique()[:5]:
        dup_data = df[df['submitter_id'] == pid]
        print(f"\n{pid} ({len(dup_data)} entries):")
        print(f"  PD-L1 values: {[f'{v:.4f}' for v in dup_data['PD1_PDL1_signature'].tolist()]}")

# 2. Clean duplicates (take mean of PD-L1)
print("\n🧼 Cleaning duplicates (Taking mean of PD1_PDL1_signature)...")

# Define aggregation functions
agg_funcs = {col: 'first' for col in df.columns if col != 'PD1_PDL1_signature'}
agg_funcs['PD1_PDL1_signature'] = 'mean'

df_clean = df.groupby('submitter_id', as_index=False).agg(agg_funcs)

print(f"✅ Cleaned data: {len(df_clean)} patients remaining (removed {dup_patients} duplicates)")

# 3. Final validation and save
print("\n🔬 Data validation (Post-cleaning):")
print(f"  Total samples: {len(df_clean)}")
print(f"  Total unique patients: {df_clean['submitter_id'].nunique()}")
print(f"  No duplicates check: {df_clean.duplicated(subset=['submitter_id']).sum() == 0}")

# Save the cleaned file
output_file = "results/immune_subtypes/immune_subtypes_TNBC_like_TCGA_core_CLEAN.tsv"
df_clean.to_csv(output_file, sep='\t', index=False)
print(f"\n💾 Saved cleaned data to: {output_file}")
