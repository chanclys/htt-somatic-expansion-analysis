import pandas as pd
import numpy as np

# Load the metrics we computed
fold_exp = pd.read_csv('HTT_fold_expansion.csv', index_col=0)
somatic_burden = pd.read_csv('HTT_somatic_burden.csv', index_col=0)
tail_metrics = pd.read_csv('HTT_tail_metrics.csv', index_col=0)

# Load raw data for additional stats
control = pd.read_csv('Control_11_HTT_detailed_metrics_v2.csv')
h1 = pd.read_csv('H1_HTT_detailed_metrics_v2.csv')
h3 = pd.read_csv('H3_HTT_detailed_metrics_v2.csv')

# Create comprehensive summary table
summary = pd.DataFrame(index=['Control_11', 'H1', 'H3'])

# Basic counts
summary['Total Reads'] = [len(control), len(h1), len(h3)]
summary['Normal Allele (AL=0) CAG'] = [
    f"{control[control['allele']==0]['CAG_count'].median():.0f} (range {control[control['allele']==0]['CAG_count'].min():.0f}-{control[control['allele']==0]['CAG_count'].max():.0f})",
    f"{h1[h1['allele']==0]['CAG_count'].median():.0f} (range {h1[h1['allele']==0]['CAG_count'].min():.0f}-{h1[h1['allele']==0]['CAG_count'].max():.0f})",
    f"{h3[h3['allele']==0]['CAG_count'].median():.0f} (range {h3[h3['allele']==0]['CAG_count'].min():.0f}-{h3[h3['allele']==0]['CAG_count'].max():.0f})",
]
summary['Expanded Allele (AL=1) CAG'] = [
    f"{control[control['allele']==1]['CAG_count'].median():.0f} (max {control[control['allele']==1]['CAG_count'].max():.0f})",
    f"{h1[h1['allele']==1]['CAG_count'].median():.0f} (max {h1[h1['allele']==1]['CAG_count'].max():.0f})",
    f"{h3[h3['allele']==1]['CAG_count'].median():.0f} (max {h3[h3['allele']==1]['CAG_count'].max():.0f})",
]

# Fold expansion
summary['Fold Expansion (median)'] = fold_exp['fold_median'].round(2)
summary['Fold Expansion (max)'] = fold_exp['fold_max'].round(1)

# Somatic burden
summary['Somatic Burden Score'] = somatic_burden['somatic_burden_score'].round(1)
summary['% Any Expansion'] = somatic_burden['pct_any_expansion'].round(1)
summary['% >100 CAG'] = somatic_burden['pct_gt_100'].round(2)
summary['% >500 CAG'] = somatic_burden['pct_gt_500'].round(2)

# Tail metrics
summary['Tail Fraction >100'] = tail_metrics['tail_fraction_gt100'].round(2)
summary['Tail Mean CAG'] = tail_metrics['tail_mean_gt100'].round(0)
summary['Skewness'] = tail_metrics['skewness'].round(2)

print("\n" + "="*150)
print("TABLE 1: HTT CAG EXPANSION METRICS - COMPREHENSIVE SUMMARY")
print("="*150)
print(summary.to_string())
print("="*150 + "\n")

summary.to_csv('HTT_SUMMARY_TABLE.csv')
print("✓ Saved HTT_SUMMARY_TABLE.csv")

# Create publication-ready text version
with open('HTT_SUMMARY_TABLE.txt', 'w') as f:
    f.write("TABLE 1: HTT CAG EXPANSION METRICS - COMPREHENSIVE SUMMARY\n")
    f.write("="*150 + "\n\n")
    f.write(summary.to_string())
    f.write("\n\n" + "="*150 + "\n")
    f.write("INTERPRETATION:\n")
    f.write("- Fold Expansion: ratio of expanded allele CAG to normal allele CAG\n")
    f.write("- Somatic Burden Score: average additional CAGs above normal baseline\n")
    f.write("- Tail Fraction >100: percentage of reads with extreme expansion (>100 CAG)\n")
    f.write("- Skewness: measure of distribution asymmetry (high skewness = long tail)\n")

print("✓ Saved HTT_SUMMARY_TABLE.txt")
