import pysam
import pandas as pd
import numpy as np

def analyze_structural_localization(bam_path, sample_name):
    """
    Characterize where structural variants are located:
    - Within HTT repeat region
    - In upstream flanking region
    - In downstream flanking region
    - Total read complexity
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
            start_offset = read.get_tag('SO')  # Start of repeat in query
            end_offset = read.get_tag('EO')    # End of repeat in query
            
            query_seq = read.query_sequence
            if not query_seq:
                continue
            
            # Define regions
            flanking_upstream_start = max(0, start_offset - 200)
            flanking_upstream_end = start_offset
            repeat_region_start = start_offset
            repeat_region_end = end_offset
            flanking_downstream_start = end_offset
            flanking_downstream_end = min(len(query_seq), end_offset + 200)
            
            # Extract sequences
            upstream_seq = query_seq[flanking_upstream_start:flanking_upstream_end] if flanking_upstream_end > flanking_upstream_start else ""
            repeat_seq = query_seq[repeat_region_start:repeat_region_end]
            downstream_seq = query_seq[flanking_downstream_start:flanking_downstream_end] if flanking_downstream_end > flanking_downstream_start else ""
            
            # Count repeats
            cag_count = repeat_seq.count('CAG') // 3
            
            # Analyze CIGAR string and map to regions
            cigar_ops = read.cigartuples if read.cigartuples else []
            
            query_pos = 0
            
            # Initialize region counters
            upstream_matches = 0
            upstream_insertions = 0
            upstream_deletions = 0
            
            repeat_matches = 0
            repeat_insertions = 0
            repeat_deletions = 0
            
            downstream_matches = 0
            downstream_insertions = 0
            downstream_deletions = 0
            
            for op_code, op_len in cigar_ops:
                # Operations that consume query sequence
                if op_code in [0, 1, 4, 7, 8]:  # M, I, S, =, X
                    op_start = query_pos
                    op_end = query_pos + op_len
                    
                    # Check overlap with each region
                    # UPSTREAM REGION
                    if op_end > flanking_upstream_start and op_start < flanking_upstream_end:
                        overlap_len = min(op_end, flanking_upstream_end) - max(op_start, flanking_upstream_start)
                        if op_code == 0:  # Match
                            upstream_matches += overlap_len
                        elif op_code == 1:  # Insertion
                            upstream_insertions += overlap_len
                        elif op_code in [4, 7, 8]:
                            upstream_matches += overlap_len
                    
                    # REPEAT REGION
                    if op_end > repeat_region_start and op_start < repeat_region_end:
                        overlap_len = min(op_end, repeat_region_end) - max(op_start, repeat_region_start)
                        if op_code == 0:  # Match
                            repeat_matches += overlap_len
                        elif op_code == 1:  # Insertion
                            repeat_insertions += overlap_len
                        elif op_code in [4, 7, 8]:
                            repeat_matches += overlap_len
                    
                    # DOWNSTREAM REGION
                    if op_end > flanking_downstream_start and op_start < flanking_downstream_end:
                        overlap_len = min(op_end, flanking_downstream_end) - max(op_start, flanking_downstream_start)
                        if op_code == 0:  # Match
                            downstream_matches += overlap_len
                        elif op_code == 1:  # Insertion
                            downstream_insertions += overlap_len
                        elif op_code in [4, 7, 8]:
                            downstream_matches += overlap_len
                    
                    query_pos += op_len
            
            # Calculate complexity scores for each region
            upstream_total = upstream_matches + upstream_insertions
            upstream_complexity = (upstream_insertions / upstream_total) if upstream_total > 0 else 0
            
            repeat_total = repeat_matches + repeat_insertions
            repeat_complexity = (repeat_insertions / repeat_total) if repeat_total > 0 else 0
            
            downstream_total = downstream_matches + downstream_insertions
            downstream_complexity = (downstream_insertions / downstream_total) if downstream_total > 0 else 0
            
            # Total read complexity
            total_insertions = upstream_insertions + repeat_insertions + downstream_insertions
            total_matches = upstream_matches + repeat_matches + downstream_matches
            total_complexity = (total_insertions / (total_matches + total_insertions)) if (total_matches + total_insertions) > 0 else 0
            
            # Classify pattern
            if repeat_insertions > upstream_insertions + downstream_insertions:
                pattern = 'Repeat-enriched'
            elif upstream_insertions + downstream_insertions > repeat_insertions * 2:
                pattern = 'Flanking-enriched'
            else:
                pattern = 'Distributed'
            
            rows.append({
                'sample': sample_name,
                'read_id': read.query_name,
                'allele': allele,
                'CAG_count': cag_count,
                # Upstream metrics
                'Upstream_insertions_bp': upstream_insertions,
                'Upstream_matches_bp': upstream_matches,
                'Upstream_complexity': upstream_complexity,
                # Repeat metrics
                'Repeat_insertions_bp': repeat_insertions,
                'Repeat_matches_bp': repeat_matches,
                'Repeat_complexity': repeat_complexity,
                # Downstream metrics
                'Downstream_insertions_bp': downstream_insertions,
                'Downstream_matches_bp': downstream_matches,
                'Downstream_complexity': downstream_complexity,
                # Total metrics
                'Total_insertions_bp': total_insertions,
                'Total_matches_bp': total_matches,
                'Total_complexity': total_complexity,
                'Structural_pattern': pattern,
                'Read_length': read.infer_query_length(),
            })
        
        except (KeyError, TypeError, AttributeError, IndexError):
            continue
    
    bamfile.close()
    return pd.DataFrame(rows)

# Analyze all samples
print("\n" + "="*100)
print("STRUCTURAL VARIANT LOCALIZATION ANALYSIS")
print("Where are the insertions? Repeat vs Flanking regions")
print("="*100 + "\n")

samples = {
    'Control_11': 'Control_11.trgt.sorted.spanning.bam',
    'H1': 'H1.trgt.sorted.spanning.bam',
    'H3': 'H3.trgt.sorted.spanning.bam',
}

all_localization_data = {}

for sample_name, bam_path in samples.items():
    print(f"Analyzing {sample_name}...")
    df = analyze_structural_localization(bam_path, sample_name)
    all_localization_data[sample_name] = df
    
    expanded = df[df['allele'] == 1]
    normal = df[df['allele'] == 0]
    
    print(f"\n  EXPANDED ALLELE STRUCTURAL BREAKDOWN:")
    print(f"    Avg CAG count: {expanded['CAG_count'].mean():.1f}")
    print(f"\n    UPSTREAM FLANKING (-200 to -1 bp from repeat):")
    print(f"      Insertions: {expanded['Upstream_insertions_bp'].mean():.2f} bp")
    print(f"      Complexity: {expanded['Upstream_complexity'].mean():.4f}")
    print(f"\n    REPEAT REGION:")
    print(f"      Insertions: {expanded['Repeat_insertions_bp'].mean():.2f} bp")
    print(f"      Complexity: {expanded['Repeat_complexity'].mean():.4f}")
    print(f"\n    DOWNSTREAM FLANKING (+1 to +200 bp from repeat):")
    print(f"      Insertions: {expanded['Downstream_insertions_bp'].mean():.2f} bp")
    print(f"      Complexity: {expanded['Downstream_complexity'].mean():.4f}")
    print(f"\n    TOTAL READ:")
    print(f"      Total insertions: {expanded['Total_insertions_bp'].mean():.2f} bp")
    print(f"      Total complexity: {expanded['Total_complexity'].mean():.4f}")
    
    print(f"\n  STRUCTURAL PATTERN DISTRIBUTION (expanded):")
    pattern_counts = expanded['Structural_pattern'].value_counts()
    for pattern, count in pattern_counts.items():
        pct = count / len(expanded) * 100
        print(f"    {pattern:20s}: {count:5d} reads ({pct:5.1f}%)")
    
    # Compare normal vs expanded
    print(f"\n  NORMAL ALLELE:")
    print(f"    Avg insertions (repeat): {normal['Repeat_insertions_bp'].mean():.2f} bp")
    print(f"    Avg insertions (total): {normal['Total_insertions_bp'].mean():.2f} bp")
    
    # Save
    df.to_csv(f'{sample_name}_structural_localization.csv', index=False)
    print(f"\n  ✓ Saved {sample_name}_structural_localization.csv\n")

# Create comprehensive comparison table
print("\n" + "="*100)
print("STRUCTURAL LOCALIZATION COMPARISON TABLE")
print("="*100 + "\n")

comparison_data = []
for sample_name, df in all_localization_data.items():
    expanded = df[df['allele'] == 1]
    
    comparison_data.append({
        'Sample': sample_name,
        'Avg_CAG': expanded['CAG_count'].mean(),
        'Upstream_ins': expanded['Upstream_insertions_bp'].mean(),
        'Repeat_ins': expanded['Repeat_insertions_bp'].mean(),
        'Downstream_ins': expanded['Downstream_insertions_bp'].mean(),
        'Total_ins': expanded['Total_insertions_bp'].mean(),
        'Repeat_%_of_total': (expanded['Repeat_insertions_bp'].sum() / expanded['Total_insertions_bp'].sum() * 100) if expanded['Total_insertions_bp'].sum() > 0 else 0,
        'Flanking_%_of_total': ((expanded['Upstream_insertions_bp'].sum() + expanded['Downstream_insertions_bp'].sum()) / expanded['Total_insertions_bp'].sum() * 100) if expanded['Total_insertions_bp'].sum() > 0 else 0,
        'Repeat-enriched_%': (expanded['Structural_pattern'] == 'Repeat-enriched').mean() * 100,
        'Flanking-enriched_%': (expanded['Structural_pattern'] == 'Flanking-enriched').mean() * 100,
        'Distributed_%': (expanded['Structural_pattern'] == 'Distributed').mean() * 100,
    })

comp_df = pd.DataFrame(comparison_data)
print(comp_df.to_string(index=False))
comp_df.to_csv('HTT_structural_localization_summary.csv', index=False)
print("\n✓ Saved HTT_structural_localization_summary.csv")

# Interpretation
print("\n\n" + "="*100)
print("INTERPRETATION GUIDE")
print("="*100 + "\n")

interpretation = """
STRUCTURAL PATTERNS:

1. REPEAT-ENRICHED:
   - Insertions concentrated within the HTT repeat itself
   - Suggests: Intramolecular duplications, repeat unit tandemizations
   - Mechanism: Slippage during replication, unequal crossing-over
   - Clinical: Intrinsic repeat instability

2. FLANKING-ENRICHED:
   - Insertions in upstream/downstream regions
   - Suggests: Transposon insertions, broad genomic rearrangements
   - Mechanism: DNA transposition, insertional mutagenesis
   - Clinical: Broader genomic instability, possibly genome-wide

3. DISTRIBUTED:
   - Insertions equally across repeat and flanking
   - Suggests: Complex rearrangements affecting large genomic segment
   - Mechanism: Multiple breakage/repair events
   - Clinical: Severe genomic instability

WHAT THIS TELLS US:
- Repeat-enriched = HTT-specific pathology
- Flanking-enriched = Systemic genomic problem (consistent with genome-wide repeat expansion you saw!)
- Distributed = Both local and systemic
"""

print(interpretation)

