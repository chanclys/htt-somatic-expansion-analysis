import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Set style for publication-quality figures
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'

def load_data(csv_files):
    """Load all detailed metrics CSVs"""
    dfs = {}
    for sample_name, csv_path in csv_files.items():
        dfs[sample_name] = pd.read_csv(csv_path)
    return dfs

def compute_fold_expansion(dfs):
    """Compute fold expansion relative to inherited normal allele"""
    results = {}
    
    for sample, df in dfs.items():
        normal_allele = df[df['allele'] == 0]
        expanded_allele = df[df['allele'] == 1]
        
        if len(normal_allele) > 0 and len(expanded_allele) > 0:
            normal_median = normal_allele['CAG_count'].median()
            expanded_median = expanded_allele['CAG_count'].median()
            expanded_max = expanded_allele['CAG_count'].max()
            
            fold_expansion_median = expanded_median / normal_median
            fold_expansion_max = expanded_max / normal_median
            
            results[sample] = {
                'normal_baseline': normal_median,
                'expanded_median': expanded_median,
                'expanded_max': expanded_max,
                'fold_median': fold_expansion_median,
                'fold_max': fold_expansion_max,
            }
    
    return pd.DataFrame(results).T

def compute_somatic_burden(dfs):
    """Compute somatic expansion burden scores"""
    results = {}
    
    for sample, df in dfs.items():
        normal_allele = df[df['allele'] == 0]
        expanded_allele = df[df['allele'] == 1]
        
        if len(normal_allele) > 0 and len(expanded_allele) > 0:
            normal_baseline = normal_allele['CAG_count'].median()
            
            # Fraction of expanded allele reads with somatic expansion
            pct_any_expansion = (expanded_allele['CAG_count'] > normal_baseline).mean() * 100
            pct_gt_50 = (expanded_allele['CAG_count'] > normal_baseline + 50).mean() * 100
            pct_gt_100 = (expanded_allele['CAG_count'] > 100).mean() * 100
            pct_gt_500 = (expanded_allele['CAG_count'] > 500).mean() * 100
            
            # Burden score: weighted by magnitude
            burden_score = (expanded_allele['CAG_count'] - normal_baseline).mean()
            
            results[sample] = {
                'pct_any_expansion': pct_any_expansion,
                'pct_gt_normal_plus_50': pct_gt_50,
                'pct_gt_100': pct_gt_100,
                'pct_gt_500': pct_gt_500,
                'somatic_burden_score': burden_score,
            }
    
    return pd.DataFrame(results).T

def compute_tail_metrics(dfs):
    """Compute tail distribution metrics (ECDF quantiles)"""
    results = {}
    
    for sample, df in dfs.items():
        expanded_allele = df[df['allele'] == 1]
        
        if len(expanded_allele) > 0:
            cag_sorted = np.sort(expanded_allele['CAG_count'].values)
            
            # Quantiles
            quantiles = [10, 25, 50, 75, 90, 95, 99]
            q_dict = {f'q{q}': np.percentile(cag_sorted, q) for q in quantiles}
            
            # Tail area (area under curve for >100 CAG)
            tail_reads = cag_sorted[cag_sorted > 100]
            tail_fraction = len(tail_reads) / len(cag_sorted) * 100
            tail_mean = tail_reads.mean() if len(tail_reads) > 0 else 0
            
            results[sample] = {
                **q_dict,
                'tail_fraction_gt100': tail_fraction,
                'tail_mean_gt100': tail_mean,
                'iqr': np.percentile(cag_sorted, 75) - np.percentile(cag_sorted, 25),
                'skewness': stats.skew(cag_sorted),
            }
    
    return pd.DataFrame(results).T

def compute_comparison_stats(dfs):
    """Compute statistical comparisons between groups"""
    print("\n" + "="*60)
    print("STATISTICAL COMPARISONS (Control vs HD samples)")
    print("="*60 + "\n")
    
    control_expanded = dfs['Control_11'][dfs['Control_11']['allele'] == 1]['CAG_count'].values
    h1_expanded = dfs['H1'][dfs['H1']['allele'] == 1]['CAG_count'].values
    h3_expanded = dfs['H3'][dfs['H3']['allele'] == 1]['CAG_count'].values
    
    # Mann-Whitney U test (non-parametric)
    u_h1, p_h1 = stats.mannwhitneyu(control_expanded, h1_expanded, alternative='less')
    u_h3, p_h3 = stats.mannwhitneyu(control_expanded, h3_expanded, alternative='less')
    u_h1h3, p_h1h3 = stats.mannwhitneyu(h1_expanded, h3_expanded, alternative='two-sided')
    
    print(f"Control_11 vs H1 (Mann-Whitney U test):")
    print(f"  U-statistic: {u_h1:.2e}")
    print(f"  p-value: {p_h1:.2e}")
    print(f"  Result: {'SIGNIFICANT' if p_h1 < 0.05 else 'NOT significant'}")
    
    print(f"\nControl_11 vs H3 (Mann-Whitney U test):")
    print(f"  U-statistic: {u_h3:.2e}")
    print(f"  p-value: {p_h3:.2e}")
    print(f"  Result: {'SIGNIFICANT' if p_h3 < 0.05 else 'NOT significant'}")
    
    print(f"\nH1 vs H3 (Mann-Whitney U test):")
    print(f"  U-statistic: {u_h1h3:.2e}")
    print(f"  p-value: {p_h1h3:.2e}")
    print(f"  Result: {'SIGNIFICANT' if p_h1h3 < 0.05 else 'NOT significant'}")

def plot_cag_distributions(dfs):
    """Plot CAG distributions by allele and sample"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    samples = ['Control_11', 'H1', 'H3']
    colors = {'allele_0': '#2E86AB', 'allele_1': '#A23B72'}
    
    for idx, sample in enumerate(samples):
        df = dfs[sample]
        
        ax = axes[idx]
        
        # Histogram for each allele
        allele_0 = df[df['allele'] == 0]['CAG_count']
        allele_1 = df[df['allele'] == 1]['CAG_count']
        
        ax.hist(allele_0, bins=30, alpha=0.6, label='Allele 0 (normal)', color=colors['allele_0'])
        ax.hist(allele_1, bins=30, alpha=0.6, label='Allele 1 (expanded)', color=colors['allele_1'])
        
        ax.set_xlabel('CAG repeat count')
        ax.set_ylabel('Number of reads')
        ax.set_title(f'{sample}\n(n={len(df)} reads)')
        ax.legend()
        ax.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig('HTT_CAG_distributions.png', dpi=300, bbox_inches='tight')
    print("✓ Saved HTT_CAG_distributions.png")
    plt.close()

def plot_ecdf(dfs):
    """Plot empirical cumulative distribution functions"""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    samples = ['Control_11', 'H1', 'H3']
    colors = ['#264653', '#2A9D8F', '#E76F51']
    
    for sample, color in zip(samples, colors):
        df = dfs[sample]
        expanded = df[df['allele'] == 1]['CAG_count'].values
        expanded_sorted = np.sort(expanded)
        
        y = np.arange(1, len(expanded_sorted) + 1) / len(expanded_sorted) * 100
        
        ax.plot(expanded_sorted, y, label=sample, linewidth=2, color=color)
    
    ax.axvline(x=100, color='red', linestyle='--', linewidth=1.5, label='100 CAG threshold', alpha=0.7)
    ax.axvline(x=500, color='orange', linestyle='--', linewidth=1.5, label='500 CAG threshold', alpha=0.7)
    
    ax.set_xlabel('CAG repeat count')
    ax.set_ylabel('Cumulative percentage (%)')
    ax.set_title('ECDF of expanded allele (AL=1) CAG counts')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('HTT_ECDF_expanded_allele.png', dpi=300, bbox_inches='tight')
    print("✓ Saved HTT_ECDF_expanded_allele.png")
    plt.close()

def plot_violin_alleles(dfs):
    """Violin plot comparing alleles across samples"""
    # Prepare data for plotting
    plot_data = []
    for sample, df in dfs.items():
        for _, row in df.iterrows():
            allele_name = f"Allele {int(row['allele'])}"
            plot_data.append({
                'Sample': sample,
                'Allele': allele_name,
                'CAG_count': row['CAG_count']
            })
    
    plot_df = pd.DataFrame(plot_data)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    sns.violinplot(data=plot_df, x='Sample', y='CAG_count', hue='Allele', ax=ax, palette=['#2E86AB', '#A23B72'])
    
    ax.set_ylabel('CAG repeat count')
    ax.set_title('CAG distributions by allele and sample')
    ax.set_ylim(0, 150)  # Focus on main range
    
    plt.tight_layout()
    plt.savefig('HTT_violin_alleles.png', dpi=300, bbox_inches='tight')
    print("✓ Saved HTT_violin_alleles.png")
    plt.close()

def plot_expansion_burden(burden_df):
    """Bar plot of somatic burden metrics"""
    fig, ax = plt.subplots(figsize=(10, 5))
    
    x = np.arange(len(burden_df))
    width = 0.25
    
    metrics = ['pct_gt_100', 'pct_gt_500', 'pct_any_expansion']
    labels = ['% reads >100 CAG', '% reads >500 CAG', '% any expansion']
    colors = ['#E76F51', '#F4A261', '#2A9D8F']
    
    for i, (metric, label, color) in enumerate(zip(metrics, labels, colors)):
        values = burden_df[metric]
        ax.bar(x + i*width, values, width, label=label, color=color)
    
    ax.set_xlabel('Sample')
    ax.set_ylabel('Percentage (%)')
    ax.set_title('Somatic expansion burden across samples')
    ax.set_xticks(x + width)
    ax.set_xticklabels(burden_df.index)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('HTT_somatic_burden.png', dpi=300, bbox_inches='tight')
    print("✓ Saved HTT_somatic_burden.png")
    plt.close()

def main():
    print("\n" + "="*60)
    print("HTT PUBLICATION-GRADE ANALYSIS")
    print("="*60)
    
    # Load data
    csv_files = {
        'Control_11': 'Control_11_HTT_detailed_metrics_v2.csv',
        'H1': 'H1_HTT_detailed_metrics_v2.csv',
        'H3': 'H3_HTT_detailed_metrics_v2.csv',
    }
    
    dfs = load_data(csv_files)
    print("\n✓ Loaded all sample data")
    
    # Compute metrics
    print("\nComputing metrics...")
    
    fold_expansion = compute_fold_expansion(dfs)
    print("\n" + "="*60)
    print("FOLD EXPANSION (relative to inherited normal allele)")
    print("="*60)
    print(fold_expansion.to_string())
    fold_expansion.to_csv('HTT_fold_expansion.csv')
    
    somatic_burden = compute_somatic_burden(dfs)
    print("\n" + "="*60)
    print("SOMATIC EXPANSION BURDEN")
    print("="*60)
    print(somatic_burden.to_string())
    somatic_burden.to_csv('HTT_somatic_burden.csv')
    
    tail_metrics = compute_tail_metrics(dfs)
    print("\n" + "="*60)
    print("TAIL DISTRIBUTION METRICS (expanded allele only)")
    print("="*60)
    print(tail_metrics.to_string())
    tail_metrics.to_csv('HTT_tail_metrics.csv')
    
    compute_comparison_stats(dfs)
    
    # Generate figures
    print("\n" + "="*60)
    print("GENERATING PUBLICATION-QUALITY FIGURES")
    print("="*60 + "\n")
    
    plot_cag_distributions(dfs)
    plot_ecdf(dfs)
    plot_violin_alleles(dfs)
    plot_expansion_burden(somatic_burden)
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE!")
    print("="*60)
    print("\nOutput files generated:")
    print("  - HTT_fold_expansion.csv")
    print("  - HTT_somatic_burden.csv")
    print("  - HTT_tail_metrics.csv")
    print("  - HTT_CAG_distributions.png")
    print("  - HTT_ECDF_expanded_allele.png")
    print("  - HTT_violin_alleles.png")
    print("  - HTT_somatic_burden.png")
    print("\n")

if __name__ == '__main__':
    main()
