import pysam
import pandas as pd
import numpy as np
from collections import Counter

def analyze_detailed_structure(bam_path, sample_name):
    """
    Detailed analysis of repeat structure and fragmentation patterns
    """
    
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    rows = []
    
    for read in bamfile:
        if not read.has_tag('TR'):
            continue
        
        try:
            repeat_id = read.get_tag('TR')
            if 'HTT' not in repeat_id:
                continue
                
            allele = read.get_tag('AL')
            start_offset = read.get_tag('SO')
            end_offset = read.get_tag('EO')
            
            query_seq = read.query_sequence
            if not query_seq:
                continue
            
            repeat_region = query_seq[max(0, start_offset):min(len(query_seq), end_offset)]
            
            # Count repeats
            cag_count = repeat_region.count('CAG') // 3
            caa_count = repeat_region.count('CAA')
            
            # Analyze CIGAR string for structural variants
            cigar_ops = read.cigartuples if read.cigartuples else []
            
            # Count operations
            num_matches = sum(op[1] for op in cigar_ops if op[0] == 0)  # M = match
            num_inserts = sum(op[1] for op in cigar_ops if op[0] == 1)  # I = insert
            num_deletes = sum(op[1] for op in cigar_ops if op[0] == 2)  # D = delete
            num_skips = sum(op[1] for op in cigar_ops if op[0] == 3)    # N = skip
            num_ops = len(cigar_ops)
            
            # Calculate complexity
            total_aligned = num_matches + num_inserts + num_deletes
            complexity_score = (num_inserts + num_deletes) / total_aligned if total_aligned > 0 else 0
            
            # Classify based on structural complexity
            if complexity_score < 0.05:
                complexity_class = 'Simple'
            elif complexity_score < 0.15:
                complexity_class = 'Moderate'
            else:
                complexity_class = 'Complex'
            
            rows.append({
                'sample': sample_name,
                'read_id': read.query_name,
                'allele': allele,
                'CAG_count': cag_count,
                'CAA_interruptions': caa_count,
                'Num_cigar_ops': num_ops,
                'Matches_bp': num_matches,
                'Insertions_bp': num_inserts,
                'Deletions_bp': num_deletes,
                'Skips_bp': num_skips,
                'Complexity_score': complexity_score,
                'Complexity_class': complexity_class,
                'Read_length': read.infer_query_length(),
                'Mapping_quality': read.mapping_quality,
                'Is_reverse': read.is_reverse,
            })
        
        except (KeyError, TypeError, AttributeError):
            continue
    
    bamfile.close()
    return pd.DataFrame(rows)

# Analyze all samples
print("\n" + "="*80)
print("DETAILED REPEAT STRUCTURE & FRAGMENTATION ANALYSIS")
print("="*80 + "\n")

samples = {
    'Control_11': 'Control_11.trgt.sorted.spanning.bam',
    'H1': 'H1.trgt.sorted.spanning.bam',
    'H3': 'H3.trgt.sorted.spanning.bam',
}

all_structure_data = {}

for sample_name, bam_path in samples.items():
    print(f"Analyzing {sample_name}...")
    df = analyze_detailed_structure(bam_path, sample_name)
    all_structure_data[sample_name] = df
    
    expanded = df[df['allele'] == 1]
    
    print(f"\n  COMPLEXITY DISTRIBUTION:")
    print(f"    Simple:    {(expanded['Complexity_class'] == 'Simple').sum():5d} reads ({(expanded['Complexity_class'] == 'Simple').mean()*100:5.1f}%)")
    print(f"    Moderate:  {(expanded['Complexity_class'] == 'Moderate').sum():5d} reads ({(expanded['Complexity_class'] == 'Moderate').mean()*100:5.1f}%)")
    print(f"    Complex:   {(expanded['Complexity_class'] == 'Complex').sum():5d} reads ({(expanded['Complexity_class'] == 'Complex').mean()*100:5.1f}%)")
    
    print(f"\n  STRUCTURAL VARIANTS (in expanded allele):")
    print(f"    Avg insertions:  {expanded['Insertions_bp'].mean():8.1f} bp")
    print(f"    Avg deletions:   {expanded['Deletions_bp'].mean():8.1f} bp")
    print(f"    Avg CIGAR ops:   {expanded['Num_cigar_ops'].mean():8.1f}")
    print(f"    Avg complexity:  {expanded['Complexity_score'].mean():8.3f}")
    
    print(f"\n  REPEAT METRICS (expanded allele):")
    print(f"    Avg CAG count:   {expanded['CAG_count'].mean():8.1f}")
    print(f"    Avg CAA ints:    {expanded['CAA_interruptions'].mean():8.1f}")
    
    # Save detailed
    df.to_csv(f'{sample_name}_detailed_structure.csv', index=False)
    print(f"\n  ✓ Saved {sample_name}_detailed_structure.csv\n")

# Compare by allele
print("\n" + "="*80)
print("STRUCTURE COMPARISON: NORMAL vs EXPANDED ALLELES")
print("="*80 + "\n")

for sample_name, df in all_structure_data.items():
    print(f"\n{sample_name}")
    print("-" * 80)
    
    normal = df[df['allele'] == 0]
    expanded = df[df['allele'] == 1]
    
    metrics = {
        'Metric': ['CAG_count', 'Insertions_bp', 'Deletions_bp', 'Complexity_score', 'Num_cigar_ops'],
        'Normal_mean': [
            normal['CAG_count'].mean(),
            normal['Insertions_bp'].mean(),
            normal['Deletions_bp'].mean(),
            normal['Complexity_score'].mean(),
            normal['Num_cigar_ops'].mean(),
        ],
        'Expanded_mean': [
            expanded['CAG_count'].mean(),
            expanded['Insertions_bp'].mean(),
            expanded['Deletions_bp'].mean(),
            expanded['Complexity_score'].mean(),
            expanded['Num_cigar_ops'].mean(),
        ]
    }
    
    comp_df = pd.DataFrame(metrics)
    comp_df['Fold_change'] = comp_df['Expanded_mean'] / comp_df['Normal_mean']
    print(comp_df.to_string(index=False))

# Create summary table
print("\n\n" + "="*80)
print("COMPLEXITY SUMMARY TABLE")
print("="*80 + "\n")

summary_data = []
for sample_name, df in all_structure_data.items():
    expanded = df[df['allele'] == 1]
    
    summary_data.append({
        'Sample': sample_name,
        'Simple_%': (expanded['Complexity_class'] == 'Simple').mean() * 100,
        'Moderate_%': (expanded['Complexity_class'] == 'Moderate').mean() * 100,
        'Complex_%': (expanded['Complexity_class'] == 'Complex').mean() * 100,
        'Avg_CAG': expanded['CAG_count'].mean(),
        'Avg_insertions': expanded['Insertions_bp'].mean(),
        'Avg_deletions': expanded['Deletions_bp'].mean(),
        'Avg_complexity_score': expanded['Complexity_score'].mean(),
    })

summary_df = pd.DataFrame(summary_data)
print(summary_df.to_string(index=False))
summary_df.to_csv('HTT_complexity_summary.csv', index=False)
print("\n✓ Saved HTT_complexity_summary.csv")

