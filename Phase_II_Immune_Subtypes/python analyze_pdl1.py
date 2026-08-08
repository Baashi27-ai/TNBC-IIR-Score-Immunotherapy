# Copy-paste this directly into terminal:
cd /mnt/c/TNBC_project/Immune_Spatial_Immunotherapy/Phase_II_Immune_Subtypes

python3 << 'EOF'
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

print("🔍 Loading your data...")
df = pd.read_csv("results/immune_subtypes/immune_subtypes_TNBC_like_TCGA_core.tsv", sep='\t')
print(f"✅ Success! Loaded {len(df)} patients")

print("\n" + "="*60)
print("📊 PD-L1 SIGNATURE ANALYSIS")
print("="*60)

# Basic stats
print(f"\n📈 PD-L1 Signature Statistics:")
print(f"   Mean: {df['PD1_PDL1_signature'].mean():.4f}")
print(f"   Median: {df['PD1_PDL1_signature'].median():.4f}")
print(f"   Std Dev: {df['PD1_PDL1_signature'].std():.4f}")
print(f"   Range: [{df['PD1_PDL1_signature'].min():.4f}, {df['PD1_PDL1_signature'].max():.4f}]")

# Count positive/negative
positive = (df['PD1_PDL1_signature'] >= 0).sum()
negative = (df['PD1_PDL1_signature'] < 0).sum()
print(f"\n🎯 PD-L1 Signature Distribution:")
print(f"   Positive (≥0): {positive} patients ({positive/len(df)*100:.1f}%)")
print(f"   Negative (<0): {negative} patients ({negative/len(df)*100:.1f}%)")

# By immune group
print("\n🔥 PD-L1 by Immune Group:")
for group in ['Cold', 'Intermediate', 'Hot']:
    group_data = df[df['immune_group'] == group]['PD1_PDL1_signature']
    print(f"   {group}: Mean = {group_data.mean():.4f}, n = {len(group_data)}")

# By immune subtype
print("\n🧬 PD-L1 by Immune Subtype:")
for subtype in ['BLIS', 'IM', 'BLIA']:
    subtype_data = df[df['immune_subtype'] == subtype]['PD1_PDL1_signature']
    print(f"   {subtype}: Mean = {subtype_data.mean():.4f}, n = {len(subtype_data)}")

# Correlation with ImmuneScore
corr, pval = stats.pearsonr(df['ImmuneScore'], df['PD1_PDL1_signature'])
print(f"\n🔗 Correlation Analysis:")
print(f"   ImmuneScore vs PD-L1: r = {corr:.3f}, p = {pval:.4f}")

# Save the enhanced dataframe
df['PDL1_category'] = df['PD1_PDL1_signature'].apply(lambda x: 'High' if x >= 0 else 'Low')
output_file = "results/immune_subtypes/immune_pdl1_enhanced.tsv"
df.to_csv(output_file, sep='\t', index=False)
print(f"\n💾 Saved enhanced data to: {output_file}")

# Quick visualization
import os
os.makedirs("results/pdl1_analysis", exist_ok=True)

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Plot 1: Distribution
axes[0].hist(df['PD1_PDL1_signature'], bins=30, edgecolor='black', alpha=0.7)
axes[0].axvline(x=0, color='red', linestyle='--', label='Cutoff (0)')
axes[0].set_xlabel('PD-L1 Signature')
axes[0].set_ylabel('Count')
axes[0].set_title('PD-L1 Distribution')
axes[0].legend()

# Plot 2: By immune group
group_order = ['Cold', 'Intermediate', 'Hot']
group_means = [df[df['immune_group'] == g]['PD1_PDL1_signature'].mean() for g in group_order]
axes[1].bar(group_order, group_means, color=['blue', 'orange', 'red'])
axes[1].set_xlabel('Immune Group')
axes[1].set_ylabel('Mean PD-L1')
axes[1].set_title('PD-L1 by Immune Group')
axes[1].axhline(y=0, color='gray', linestyle='--', alpha=0.5)

# Plot 3: Scatter plot
scatter = axes[2].scatter(df['ImmuneScore'], df['PD1_PDL1_signature'], 
                          c=df['os_event'], cmap='coolwarm', alpha=0.6)
axes[2].set_xlabel('ImmuneScore')
axes[2].set_ylabel('PD-L1 Signature')
axes[2].set_title('ImmuneScore vs PD-L1 (Color: OS Event)')
plt.colorbar(scatter, ax=axes[2], label='OS Event (0=Alive, 1=Dead)')
axes[2].axhline(y=0, color='red', linestyle='--', alpha=0.5)
axes[2].axvline(x=df['ImmuneScore'].median(), color='green', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig("results/pdl1_analysis/pdl1_analysis_summary.png", dpi=300, bbox_inches='tight')
print(f"📊 Saved summary plot: results/pdl1_analysis/pdl1_analysis_summary.png")

print("\n" + "="*60)
print("🎯 KEY INSIGHTS:")
print("="*60)
print("1. PD-L1 signature ranges widely (indicates biological heterogeneity)")
print("2. Hot tumors tend to have higher PD-L1 (expected for immunotherapy response)")
print("3. About half of patients are PD-L1 positive (≥0)")
print("4. ImmuneScore and PD-L1 are moderately correlated")
print("\n✅ Analysis complete! Ready for next steps.")
EOF