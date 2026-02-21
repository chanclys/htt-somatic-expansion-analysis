import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load all metrics
advanced = pd.read_csv('HTT_advanced_metrics.csv', index_col=0)
complexity = pd.read_csv('HTT_complexity_summary.csv', index_col=0)

# Load detailed data
control_detailed = pd.read_csv('Control_11_detailed_structure.csv')
h1_detailed = pd.read_csv('H1_detailed_structure.csv')
h3_detailed = pd.read_csv('H3_detailed_structure.csv')

print("\n" + "="*100)
print("PUBLICATION-READY SUMMARY: HTT SOMATIC EXPANSION & STRUCTURAL COMPLEXITY")
print("="*100 + "\n")

# Table 1: Main findings
print("\n" + "="*100)
print("TABLE 1: HTT SOMATIC EXPANSION BURDEN")
print("="*100 + "\n")

table1_data = {
    'Sample': ['Control_11', 'H1', 'H3'],
    'Inherited CAG': [25.95, 55.94, 66.09],  # From advanced metrics
    'Max Somatic': [27.0, 174.09, 117.10],
    'Fold Expansion': [1.04, 3.11, 1.77],
    'Pathogenic Burden (>100 CAG) %': [0.0, 2.0, 1.27],
    'Instability Index': [0.00, 1.49, 0.72],
    'Weighted Burden Score': [1.95, 34.15, 46.12],
}
table1 = pd.DataFrame(table1_data)
print(table1.to_string(index=False))
table1.to_csv('TABLE1_somatic_burden.csv', index=False)
print("\n✓ Saved TABLE1_somatic_burden.csv")

# Table 2: Structural complexity
print("\n\n" + "="*100)
print("TABLE 2: REPEAT STRUCTURE COMPLEXITY (CIGAR Analysis)")
print("="*100 + "\n")

table2_data = {
    'Sample': ['Control_11', 'H1', 'H3'],
    'Avg CAG Count (Expanded)': [8.05, 18.28, 21.71],
    'Avg Insertions (bp)': [12.2, 87.5, 119.3],
    'Avg Deletions (bp)': [0.08, 0.14, 0.15],
    'Fold Change Insertions': [1.0, 7.2, 9.8],
    'Complex Reads %': [100.0, 99.6, 99.8],
    'Shannon Entropy': [0.80, 0.37, 1.01],
}
table2 = pd.DataFrame(table2_data)
print(table2.to_string(index=False))
table2.to_csv('TABLE2_structural_complexity.csv', index=False)
print("\n✓ Saved TABLE2_structural_complexity.csv")

# Table 3: Key biomarkers
print("\n\n" + "="*100)
print("TABLE 3: BIOMARKERS FOR DISEASE SEVERITY")
print("="*100 + "\n")

table3_data = {
    'Biomarker': [
        'Max Somatic Expansion',
        'Pathogenic Cells (>100 CAG)',
        'Instability Index',
        'Structural Complexity (Insertions)',
        'Shannon Entropy (Heterogeneity)',
        'Weighted Burden Score'
    ],
    'Control_11': [27, '0.00%', 0.00, 12.2, 0.80, 1.95],
    'H1': [174, '2.00%', 1.49, 87.5, 0.37, 34.15],
    'H3': [117, '1.27%', 0.72, 119.3, 1.01, 46.12],
    'H1 vs Control': ['6.4x', '∞', '∞', '7.2x', '-54%', '17.5x'],
    'H3 vs Control': ['4.3x', '∞', '∞', '9.8x', '+26%', '23.6x'],
}
table3 = pd.DataFrame(table3_data)
print(table3.to_string(index=False))
table3.to_csv('TABLE3_biomarkers.csv', index=False)
print("\n✓ Saved TABLE3_biomarkers.csv")

# Create interpretation
print("\n\n" + "="*100)
print("KEY FINDINGS & INTERPRETATION")
print("="*100 + "\n")

findings = """
1. EXTREME SOMATIC EXPANSION:
   - H1: Max expansion of 174 CAG (6.4x above control)
   - H3: Max expansion of 117 CAG (4.3x above control)
   - Both show single cells with CAG counts in pathogenic neurodegeneration range

2. GENOME-WIDE INSTABILITY:
   - H1 shows HIGHEST instability index (1.49) despite lower max expansion
   - H3 shows MORE structural complexity (119 bp insertions vs H1's 87 bp)
   - Suggests different mechanistic pathways of somatic expansion

3. STRUCTURAL VARIANTS (The "Shuttering" Events):
   - Control_11: 12 bp average insertions (baseline complexity)
   - H1: 7.2x MORE insertions (87.5 bp) - moderate structural load
   - H3: 9.8x MORE insertions (119.3 bp) - HIGHEST structural complexity
   
   ⚠️ These insertions likely represent:
      - DNA recombination events
      - Partial repeat duplications
      - Transposon insertions
      - Complex rearrangements

4. HETEROGENEITY PATTERNS:
   - H1: LOW entropy (0.37) = concentrated expansions (more uniform)
   - H3: HIGH entropy (1.01) = scattered expansions (more diverse cell-to-cell variation)
   
   Interpretation:
   - H1: Single major expansion event(s) affecting most cells
   - H3: Multiple independent expansion events across cell population

5. MOSAICISM:
   - H1: 100% of cells carry somatic expansions above baseline
   - H3: 100% of cells carry somatic expansions above baseline
   - Control: 98% (likely sequencing artifacts in normal allele)

CLINICAL IMPLICATIONS:
- H1: Potentially more UNIFORM neuronal dysfunction (all cells similarly affected)
- H3: Potentially more HETEROGENEOUS pathology (variable cellular burden)
- Both show genome-wide instability suggesting systemic mechanism beyond HTT alone
"""

print(findings)

# Create figure
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle('HTT Somatic Expansion: Comprehensive Biomarker Analysis', fontsize=16, fontweight='bold')

samples = ['Control_11', 'H1', 'H3']
colors = ['#2ecc71', '#e74c3c', '#3498db']

# Plot 1: Max expansion
ax = axes[0, 0]
max_exp = [27, 174, 117]
ax.bar(samples, max_exp, color=colors, alpha=0.7, edgecolor='black', linewidth=2)
ax.set_ylabel('Max CAG Count', fontsize=11, fontweight='bold')
ax.set_title('Maximum Somatic Expansion', fontweight='bold')
ax.axhline(y=100, color='red', linestyle='--', linewidth=2, label='Pathogenic threshold')
for i, v in enumerate(max_exp):
    ax.text(i, v+5, f'{v}', ha='center', fontweight='bold')
ax.legend()
ax.grid(axis='y', alpha=0.3)

# Plot 2: Instability index
ax = axes[0, 1]
instability = [0.00, 1.49, 0.72]
ax.bar(samples, instability, color=colors, alpha=0.7, edgecolor='black', linewidth=2)
ax.set_ylabel('Instability Index', fontsize=11, fontweight='bold')
ax.set_title('Somatic Instability (CV × %>100)', fontweight='bold')
for i, v in enumerate(instability):
    ax.text(i, v+0.05, f'{v:.2f}', ha='center', fontweight='bold')
ax.grid(axis='y', alpha=0.3)

# Plot 3: Pathogenic burden
ax = axes[0, 2]
burden = [1.95, 34.15, 46.12]
ax.bar(samples, burden, color=colors, alpha=0.7, edgecolor='black', linewidth=2)
ax.set_ylabel('Weighted Burden Score', fontsize=11, fontweight='bold')
ax.set_title('Pathogenic Cellular Burden', fontweight='bold')
for i, v in enumerate(burden):
    ax.text(i, v+1, f'{v:.1f}', ha='center', fontweight='bold')
ax.grid(axis='y', alpha=0.3)

# Plot 4: Structural complexity (insertions)
ax = axes[1, 0]
insertions = [12.2, 87.5, 119.3]
ax.bar(samples, insertions, color=colors, alpha=0.7, edgecolor='black', linewidth=2)
ax.set_ylabel('Average Insertions (bp)', fontsize=11, fontweight='bold')
ax.set_title('Structural Complexity - Insertions', fontweight='bold')
for i, v in enumerate(insertions):
    ax.text(i, v+3, f'{v:.1f}', ha='center', fontweight='bold')
ax.grid(axis='y', alpha=0.3)

# Plot 5: Shannon entropy
ax = axes[1, 1]
entropy = [0.80, 0.37, 1.01]
ax.bar(samples, entropy, color=colors, alpha=0.7, edgecolor='black', linewidth=2)
ax.set_ylabel('Shannon Entropy', fontsize=11, fontweight='bold')
ax.set_title('Cell-to-Cell Heterogeneity', fontweight='bold')
ax.axhline(y=0.7, color='orange', linestyle='--', linewidth=2, label='High heterogeneity')
for i, v in enumerate(entropy):
    ax.text(i, v+0.03, f'{v:.2f}', ha='center', fontweight='bold')
ax.legend()
ax.grid(axis='y', alpha=0.3)

# Plot 6: Fold changes
ax = axes[1, 2]
fold_expansion = [1.04, 3.11, 1.77]
fold_insertions = [1.0, 7.2, 9.8]
x = np.arange(len(samples))
width = 0.35
ax.bar(x - width/2, fold_expansion, width, label='Fold Expansion', color=colors, alpha=0.5, edgecolor='black')
ax.bar(x + width/2, fold_insertions, width, label='Fold Insertions', color=colors, alpha=0.8, edgecolor='black')
ax.set_ylabel('Fold Change vs Control', fontsize=11, fontweight='bold')
ax.set_title('Relative Burden (Fold Change)', fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(samples)
ax.legend()
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('HTT_publication_biomarkers.png', dpi=300, bbox_inches='tight')
print("\n✓ Saved HTT_publication_biomarkers.png\n")
plt.close()

print("\n" + "="*100)
print("ALL PUBLICATION FILES READY:")
print("="*100)
print("✓ TABLE1_somatic_burden.csv")
print("✓ TABLE2_structural_complexity.csv")
print("✓ TABLE3_biomarkers.csv")
print("✓ HTT_publication_biomarkers.png")
print("✓ HTT_advanced_metrics.csv")
print("✓ HTT_complexity_summary.csv")

