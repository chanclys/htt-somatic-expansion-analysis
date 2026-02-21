import pandas as pd
import numpy as np
from scipy import stats

# Load the genome-wide data
gw_summary = pd.read_csv('genome_wide_repeats_summary.csv')

print("\n" + "="*80)
print("IDENTIFYING ABNORMAL REPEAT LOCI (Potential Genome-Wide Instability)")
print("="*80 + "\n")

# Pivot to compare samples
pivot_df = gw_summary.pivot(index='Repeat_ID', columns='Sample', values='Avg_repeat_units')

# Calculate differences
pivot_df['H1_vs_Control'] = pivot_df['H1'] - pivot_df['Control_11']
pivot_df['H3_vs_Control'] = pivot_df['H3'] - pivot_df['Control_11']
pivot_df['H1_vs_H3'] = pivot_df['H1'] - pivot_df['H3']

# Calculate fold changes
pivot_df['H1_fold_change'] = pivot_df['H1'] / pivot_df['Control_11']
pivot_df['H3_fold_change'] = pivot_df['H3'] / pivot_df['Control_11']

# Mark repeats with significant expansions
pivot_df['H1_expanded'] = (pivot_df['H1_fold_change'] > 1.2) | (pivot_df['H1_vs_Control'] > 10)
pivot_df['H3_expanded'] = (pivot_df['H3_fold_change'] > 1.2) | (pivot_df['H3_vs_Control'] > 10)

# Sort by H1 fold change
pivot_df_sorted = pivot_df.sort_values('H1_fold_change', ascending=False)

print("TOP 15 REPEATS WITH LARGEST EXPANSIONS IN H1 vs CONTROL:")
print("-" * 80)
top_h1 = pivot_df_sorted[['Control_11', 'H1', 'H3', 'H1_fold_change', 'H1_vs_Control']].head(15)
print(top_h1.to_string())

print("\n\nTOP 15 REPEATS WITH LARGEST EXPANSIONS IN H3 vs CONTROL:")
print("-" * 80)
pivot_df_h3 = pivot_df.sort_values('H3_fold_change', ascending=False)
top_h3 = pivot_df_h3[['Control_11', 'H1', 'H3', 'H3_fold_change', 'H3_vs_Control']].head(15)
print(top_h3.to_string())

# Identify repeats that are abnormal in BOTH HD patients
print("\n\nREPEATS ABNORMAL IN BOTH H1 AND H3 (potential biomarkers):")
print("-" * 80)
both_expanded = pivot_df[(pivot_df['H1_expanded']) & (pivot_df['H3_expanded'])].copy()
both_expanded_sorted = both_expanded.sort_values('H1_fold_change', ascending=False)
print(both_expanded_sorted[['Control_11', 'H1', 'H3', 'H1_fold_change', 'H3_fold_change']].to_string())

# Statistical summary
print("\n\n" + "="*80)
print("SUMMARY STATISTICS")
print("="*80)
print(f"\nTotal repeat loci analyzed: {len(pivot_df)}")
print(f"Repeats expanded in H1: {pivot_df['H1_expanded'].sum()}")
print(f"Repeats expanded in H3: {pivot_df['H3_expanded'].sum()}")
print(f"Repeats expanded in BOTH: {(pivot_df['H1_expanded'] & pivot_df['H3_expanded']).sum()}")

# Save detailed results
pivot_df.to_csv('HTT_genome_wide_repeat_comparison.csv')
print(f"\n✓ Saved HTT_genome_wide_repeat_comparison.csv")

# Create focused table on HTT and nearby repeats
print("\n" + "="*80)
print("HTT AND RELATED LOCI")
print("="*80 + "\n")
htt_related = ['HD_HTT', 'HAP1', 'FMR1', 'FXS_FMR1', 'FRDA_FXN', 'DM1_DMPK']
htt_table = pivot_df.loc[pivot_df.index.isin(htt_related), ['Control_11', 'H1', 'H3', 'H1_fold_change', 'H3_fold_change']]
if len(htt_table) > 0:
    print(htt_table.to_string())
    htt_table.to_csv('HTT_related_loci.csv')
    print(f"\n✓ Saved HTT_related_loci.csv")

