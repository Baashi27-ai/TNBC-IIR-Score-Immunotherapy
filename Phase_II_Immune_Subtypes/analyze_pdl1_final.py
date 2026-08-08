import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import warnings
import os
warnings.filterwarnings('ignore') # Suppress warnings for cleaner output

# ==================== CONFIGURATION ====================
# Load the CLEANED file (Corrected input path)
INPUT_FILE = "results/immune_subtypes/immune_subtypes_TNBC_like_TCGA_core_CLEAN.tsv"
OUTPUT_DIR = "results/pdl1_analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ==================== LOAD DATA ====================
print("📊 Loading CLEANED data for final analysis...")
try:
    df = pd.read_csv(INPUT_FILE, sep='\t')
    print(f"✅ Clean data loaded: {len(df)} patients")
except FileNotFoundError:
    print(f"❌ Error: Cleaned input file not found at: {INPUT_FILE}")
    print("Please ensure the data cleaning step was run successfully.")
    exit(1)


# ==================== QUICK INSIGHTS & INITIAL GROUPING ====================
print("\n" + "="*60)
print("📊 QUICK INSIGHTS FROM CLEAN DATA")
print("="*60)

# 1. Basic PD-L1 stats
print(f"\n📈 PD-L1 Signature:")
print(f"   Mean: {df['PD1_PDL1_signature'].mean():.3f}")
print(f"   Median: {df['PD1_PDL1_signature'].median():.3f}")
print(f"   Range: [{df['PD1_PDL1_signature'].min():.3f}, {df['PD1_PDL1_signature'].max():.3f}]")

# 2. Categorize PD-L1 (High/Low)
pdl1_median = df['PD1_PDL1_signature'].median()
df['PDL1_group'] = df['PD1_PDL1_signature'].apply(
    lambda x: 'High' if x >= pdl1_median else 'Low'
)

print(f"\n🎯 PD-L1 Groups (median cutoff = {pdl1_median:.3f}):")
print(f"   High: {(df['PDL1_group'] == 'High').sum()} patients ({(df['PDL1_group'] == 'High').sum()/len(df)*100:.1f}%)")
print(f"   Low: {(df['PDL1_group'] == 'Low').sum()} patients ({(df['PDL1_group'] == 'Low').sum()/len(df)*100:.1f}%)")

# 3. PD-L1 by Immune Subtype
print("\n🔥 PD-L1 by Immune Subtype:")
for subtype in ['BLIA', 'IM', 'BLIS']:
    subtype_data = df[df['immune_subtype'] == subtype]
    mean_pdl1 = subtype_data['PD1_PDL1_signature'].mean()
    high_pct = (subtype_data['PDL1_group'] == 'High').mean() * 100
    print(f"   {subtype}: Mean PD-L1 = {mean_pdl1:.3f}, {high_pct:.1f}% High")

# 4. Correlation with ImmuneScore
corr, pval = stats.pearsonr(df['ImmuneScore'], df['PD1_PDL1_signature'])
print(f"\n🔗 Correlation: ImmuneScore vs PD-L1 = {corr:.3f} (p = {pval:.4f})")

# 5. Survival stats
print(f"\n💀 Survival Events: {df['os_event'].sum()} deaths out of {len(df)} patients ({(df['os_event'].sum()/len(df))*100:.1f}%)")

# ==================== VISUALIZATION (Comprehensive Plot) ====================
print("\n🎨 Creating comprehensive visualizations...")
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Plot 1: PD-L1 distribution
axes[0,0].hist(df['PD1_PDL1_signature'], bins=30, edgecolor='black', alpha=0.7, color='skyblue')
axes[0,0].axvline(x=pdl1_median, color='red', linestyle='--', label=f'Median = {pdl1_median:.2f}')
axes[0,0].set_xlabel('PD-L1 Signature')
axes[0,0].set_ylabel('Count')
axes[0,0].set_title('A) PD-L1 Distribution')
axes[0,0].legend()
axes[0,0].grid(True, alpha=0.3)

# Plot 2: PD-L1 by immune subtype
sns.boxplot(data=df, x='immune_subtype', y='PD1_PDL1_signature', 
             order=['BLIA', 'IM', 'BLIS'], ax=axes[0,1])
axes[0,1].axhline(y=pdl1_median, color='red', linestyle='--', alpha=0.5)
axes[0,1].set_xlabel('Immune Subtype')
axes[0,1].set_ylabel('PD-L1 Signature')
axes[0,1].set_title('B) PD-L1 by Immune Subtype')
axes[0,1].grid(True, alpha=0.3)


# Plot 3: ImmuneScore vs PD-L1
scatter = axes[0,2].scatter(df['ImmuneScore'], df['PD1_PDL1_signature'], 
                            c=df['os_event'], cmap='coolwarm', alpha=0.6, s=50)
axes[0,2].set_xlabel('ImmuneScore')
axes[0,2].set_ylabel('PD-L1 Signature')
axes[0,2].set_title(f'C) ImmuneScore vs PD-L1 (r = {corr:.3f})')
plt.colorbar(scatter, ax=axes[0,2], label='OS Event (0=Alive, 1=Dead)')
axes[0,2].axhline(y=pdl1_median, color='red', linestyle='--', alpha=0.5)
axes[0,2].grid(True, alpha=0.3)

# Plot 4: Placeholder for KM curve (will be plotted separately)
axes[1,0].set_axis_off() 

# Plot 5: Stage distribution by PD-L1
# Robustly clean stage names for plotting
df['stage_clean'] = df['ajcc_pathologic_stage'].fillna('N/A').astype(str).str.replace('Stage ', '').str.strip()
stage_counts = pd.crosstab(df['stage_clean'], df['PDL1_group'])
# Define order for better visualization
stage_order = ['IA', 'IB', 'I', 'IIA', 'IIB', 'II', 'IIIA', 'IIIB', 'IIIC', 'III', 'IV', 'N/A']
stage_counts = stage_counts.reindex(stage_order, fill_value=0).iloc[::] 

stage_counts.plot(kind='bar', ax=axes[1,1])
axes[1,1].set_xlabel('AJCC Pathologic Stage')
axes[1,1].set_ylabel('Count')
axes[1,1].set_title('E) Stage Distribution by PD-L1 Group')
axes[1,1].legend(title='PD-L1 Group')
axes[1,1].tick_params(axis='x', rotation=45)
axes[1,1].grid(True, alpha=0.3)

# Plot 6: Combined Immuno-Score
# Create a combined score (Re-run normalization on cleaned data)
df['ImmunoScore_norm'] = (df['ImmuneScore'] - df['ImmuneScore'].min()) / (df['ImmuneScore'].max() - df['ImmuneScore'].min())
df['PDL1_norm'] = (df['PD1_PDL1_signature'] - df['PD1_PDL1_signature'].min()) / (df['PD1_PDL1_signature'].max() - df['PD1_PDL1_signature'].min())
df['Combined_ImmunoScore'] = 0.6 * df['ImmunoScore_norm'] + 0.4 * df['PDL1_norm']

axes[1,2].hist(df['Combined_ImmunoScore'], bins=30, edgecolor='black', alpha=0.7, color='green')
axes[1,2].set_xlabel('Combined ImmunoScore')
axes[1,2].set_ylabel('Count')
axes[1,2].set_title('F) Combined Immune + PD-L1 Score')
axes[1,2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/pdl1_comprehensive_analysis.png", dpi=300, bbox_inches='tight')
print(f"✅ Saved: {OUTPUT_DIR}/pdl1_comprehensive_analysis.png")


# ==================== SURVIVAL ANALYSIS (KM Plots) ====================
print("\n📊 Running survival analysis...")

# Remove impossible OS times (negative, too short)
df_survival = df[df['os_time'].notna() & (df['os_time'] > 30)].copy()
print(f"   Using {len(df_survival)} patients for survival analysis")

plt.figure(figsize=(16, 8))

# Plot 1 (Matches placeholder D): PD-L1 High vs Low
plt.subplot(1, 2, 1)
ax1 = plt.gca()
colors = {'High': 'red', 'Low': 'blue'}
for group in ['High', 'Low']:
    mask = df_survival['PDL1_group'] == group
    kmf = KaplanMeierFitter()
    kmf.fit(df_survival.loc[mask, 'os_time'], 
             df_survival.loc[mask, 'os_event'], 
             label=f'PD-L1 {group} (n={mask.sum()})')
    kmf.plot_survival_function(ax=ax1, color=colors[group])

# Log-rank test
results = logrank_test(
    df_survival.loc[df_survival['PDL1_group'] == 'High', 'os_time'],
    df_survival.loc[df_survival['PDL1_group'] == 'Low', 'os_time'],
    df_survival.loc[df_survival['PDL1_group'] == 'High', 'os_event'],
    df_survival.loc[df_survival['PDL1_group'] == 'Low', 'os_event']
)

plt.title(f'D) Survival by PD-L1 Signature (Log-rank p = {results.p_value:.4f})')
plt.xlabel('Time (days)')
plt.ylabel('Survival Probability')
plt.grid(True, alpha=0.3)
plt.legend()


# Plot 2: Combined groups (Immunosubtype + PD-L1)
plt.subplot(1, 2, 2)
# Create combined groups: Hot/High, Hot/Low, Cold/High, Cold/Low
df_survival['combined_group'] = df_survival['immune_group'] + '/' + df_survival['PDL1_group']
group_order = ['Hot/High', 'Hot/Low', 'Intermediate/High', 'Intermediate/Low', 'Cold/High', 'Cold/Low']

colors_combined = {'Hot/High': 'darkred', 'Hot/Low': 'lightcoral',
                   'Intermediate/High': 'darkorange', 'Intermediate/Low': 'navajowhite',
                   'Cold/High': 'darkblue', 'Cold/Low': 'lightblue'}

for group in group_order:
    mask = df_survival['combined_group'] == group
    if mask.sum() > 0:  # Only plot if group exists
        kmf = KaplanMeierFitter()
        kmf.fit(df_survival.loc[mask, 'os_time'], 
                df_survival.loc[mask, 'os_event'], 
                label=f'{group} (n={mask.sum()})')
        kmf.plot_survival_function(ax=plt.gca(), color=colors_combined.get(group, 'gray'))

plt.title('G) Survival by Combined Immune+PD-L1 Groups')
plt.xlabel('Time (days)')
plt.ylabel('Survival Probability')
plt.grid(True, alpha=0.3)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/survival_analysis.png", dpi=300, bbox_inches='tight')
print(f"✅ Saved: {OUTPUT_DIR}/survival_analysis.png")


# ==================== SAVE FINAL ENHANCED DATA & SUMMARY ====================
print("\n💾 Saving enhanced data...")

# Add therapy recommendations
def recommend_therapy(row):
    if row['immune_subtype'] == 'BLIA' and row['PDL1_group'] == 'High':
        return "STRONG candidate for anti-PD1/L1 therapy"
    elif row['immune_subtype'] == 'BLIA' and row['PDL1_group'] == 'Low':
        return "Consider anti-PD1/L1 + chemotherapy combo"
    elif row['immune_subtype'] == 'BLIS' and row['PDL1_group'] == 'High':
        return "Consider clinical trials: TME modulation + immunotherapy"
    elif row['immune_subtype'] == 'BLIS' and row['PDL1_group'] == 'Low':
        return "Standard chemotherapy ± novel combinations"
    elif row['immune_subtype'] == 'IM':
        return "Consider biomarker-guided combination therapy"
    else:
        return "Individualized therapy consideration needed"

df['Therapy_Recommendation'] = df.apply(recommend_therapy, axis=1)

# Save enhanced data
output_file = f"{OUTPUT_DIR}/immune_pdl1_enhanced.tsv"
df.to_csv(output_file, sep='\t', index=False)
print(f"✅ Saved enhanced data: {output_file}")

# Summary statistics
print("\n" + "="*60)
print("🎯 CLINICAL IMPLICATIONS SUMMARY")
print("="*60)

# Count therapy recommendations
therapy_counts = df['Therapy_Recommendation'].value_counts()
print("\n💊 Therapy Recommendations:")
for therapy, count in therapy_counts.items():
    percent = count/len(df)*100
    print(f"  {therapy}: {count} patients ({percent:.1f}%)")

# Best candidates for immunotherapy
best_candidates = df[(df['immune_subtype'] == 'BLIA') & (df['PDL1_group'] == 'High')]
print(f"\n🔥 BEST immunotherapy candidates (BLIA + High PD-L1):")
print(f"   {len(best_candidates)} patients ({len(best_candidates)/len(df)*100:.1f}% of cohort)")

# Worst candidates
worst_candidates = df[(df['immune_subtype'] == 'BLIS') & (df['PDL1_group'] == 'Low')]
print(f"\n❌ POOR immunotherapy candidates (BLIS + Low PD-L1):")
print(f"   {len(worst_candidates)} patients ({len(worst_candidates)/len(df)*100:.1f}% of cohort)")

print("\n✅ Final Analysis complete!")
print(f"📁 Results saved in: {OUTPUT_DIR}/")
