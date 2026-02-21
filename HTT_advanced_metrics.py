import pandas as pd
import numpy as np
from scipy import stats

# Load data
control = pd.read_csv('Control_11_HTT_detailed_metrics_v2.csv')
h1 = pd.read_csv('H1_HTT_detailed_metrics_v2.csv')
h3 = pd.read_csv('H3_HTT_detailed_metrics_v2.csv')

print("\n" + "="*80)
print("ADVANCED METRICS FOR PUBLICATION")
print("="*80 + "\n")

def compute_advanced_metrics(df, sample_name):
    """Compute advanced somatic expansion metrics"""
    
    expanded = df[df['allele'] == 1]['CAG_count'].values
    normal = df[df['allele'] == 0]['CAG_count'].values
    
    normal_baseline = np.median(normal)
    
    metrics = {}
    
    # 1. DISTRIBUTION SHAPE METRICS
    metrics['Mean_expanded'] = expanded.mean()
    metrics['Median_expanded'] = np.median(expanded)
    metrics['Std_expanded'] = expanded.std()
    metrics['CV_expanded'] = expanded.std() / expanded.mean()  # Coefficient of variation
    metrics['Range_expanded'] = expanded.max() - expanded.min()
    
    # 2. SOMATIC EXPANSION SEVERITY
    metrics['Max_expansion'] = expanded.max()
    metrics['Mean_expansion_above_baseline'] = (expanded - normal_baseline).mean()
    metrics['Median_expansion_above_baseline'] = np.median(expanded - normal_baseline)
    
    # 3. TAIL METRICS (pathogenic range)
    metrics['Pct_gt_100'] = (expanded > 100).mean() * 100
    metrics['Pct_gt_200'] = (expanded > 200).mean() * 100
    metrics['Pct_gt_500'] = (expanded > 500).mean() * 100
    metrics['Pct_gt_1000'] = (expanded > 1000).mean() * 100
    
    # 4. PERCENTILE-BASED METRICS
    percentiles = [5, 10, 25, 50, 75, 90, 95, 99]
    for p in percentiles:
        metrics[f'P{p}'] = np.percentile(expanded, p)
    
    # 5. DISPERSION METRICS
    q1 = np.percentile(expanded, 25)
    q3 = np.percentile(expanded, 75)
    metrics['IQR'] = q3 - q1
    metrics['QCD'] = (q3 - q1) / (q3 + q1) if (q3 + q1) > 0 else 0  # Quartile coefficient of dispersion
    
    # 6. ENTROPY/DIVERSITY
    hist, _ = np.histogram(expanded, bins=50)
    hist = hist[hist > 0]  # Remove zeros
    probs = hist / hist.sum()
    metrics['Shannon_entropy'] = -np.sum(probs * np.log2(probs))  # Measure of distribution diversity
    
    # 7. SOMATIC BURDEN INDEX
    # Weighted by severity: reads >100 get higher weight
    burden_weights = np.where(expanded > 100, 2.0, 1.0)
    metrics['Weighted_burden_score'] = ((expanded - normal_baseline) * burden_weights).mean()
    
    # 8. INSTABILITY INDEX
    # Higher = more unstable (more variation in expansion)
    metrics['Instability_index'] = metrics['CV_expanded'] * metrics['Pct_gt_100']
    
    # 9. MOSAICISM RATIO
    # Ratio of cells with >inherited baseline expansion
    metrics['Mosaicism_ratio'] = (expanded > normal_baseline).sum() / len(expanded)
    
    # 10. PATHOGENIC LOAD
    # Sum of all expansions above pathogenic threshold (100 CAG)
    pathogenic_threshold = 100
    pathogenic_reads = expanded[expanded > pathogenic_threshold]
    metrics['Pathogenic_load'] = (pathogenic_reads - pathogenic_threshold).sum() if len(pathogenic_reads) > 0 else 0
    
    return metrics

# Compute for all samples
samples = {'Control_11': control, 'H1': h1, 'H3': h3}
all_metrics = {}

for sample_name, df in samples.items():
    metrics = compute_advanced_metrics(df, sample_name)
    all_metrics[sample_name] = metrics
    
    print(f"\n{sample_name}")
    print("-" * 80)
    for key, value in sorted(metrics.items()):
        if isinstance(value, float):
            print(f"  {key:40s}: {value:15.2f}")

# Create comparison table
metrics_df = pd.DataFrame(all_metrics).T
metrics_df.to_csv('HTT_advanced_metrics.csv')
print("\n✓ Saved HTT_advanced_metrics.csv")

# Create interpretation guide
print("\n" + "="*80)
print("METRIC INTERPRETATION GUIDE")
print("="*80 + "\n")

interpretation = """
DISTRIBUTION SHAPE:
  - Mean/Median: Central tendency of expansion
  - Std/CV: Variability (higher = more heterogeneous)
  - Range: Span from smallest to largest expansion

SEVERITY METRICS:
  - Max_expansion: Largest observed somatic expansion
  - Mean/Median_expansion_above_baseline: Average somatic burden

TAIL METRICS (% reads with extreme expansion):
  - Pct_gt_100/200/500/1000: Fraction of cells with pathogenic expansions
  - Higher % = more cells at risk for neurodegeneration

PERCENTILES:
  - P50 = Median (50th percentile)
  - P90/P95/P99 = Extreme expansions (top 10%, 5%, 1% of cells)

DISPERSION:
  - IQR: Interquartile range (middle 50% spread)
  - QCD: Normalized dispersion (0-1 scale)

SHANNON ENTROPY:
  - Measure of distribution diversity
  - Higher = more "spread out" (more cell-to-cell variation)
  - Lower = more "concentrated" (more uniform expansion)

SOMATIC BURDEN INDEX:
  - Weighted by severity (extreme expansions count more)
  - Better capture of pathogenic load than simple average

INSTABILITY INDEX:
  - Combines coefficient of variation with % extreme expansions
  - Higher = more unstable (unpredictable somatic expansion)
  - Useful biomarker for disease progression

MOSAICISM RATIO:
  - Fraction of cells with ANY expansion above inherited baseline
  - 1.0 = 100% of cells are somatic mosaics
  - Important for understanding tissue burden

PATHOGENIC LOAD:
  - Sum of all expansions above 100 CAG threshold
  - Captures total "excess" CAG repeats in pathogenic range
  - Correlates with cellular stress
"""

print(interpretation)

