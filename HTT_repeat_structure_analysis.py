import pysam
import pandas as pd
import numpy as np

def analyze_repeat_structure(bam_path, sample_name):
    """
    Analyze repeat structure variants:
    - Clean CAG repeats
    - CAA interruptions
    - Fragmented repeats (gaps in reads)
    - Repeat unit composition
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
            
            # CLASSIFY REPEAT STRUCTURE
            cag_count = repeat_region.count('CAG') // 3
            caa_count = repeat_region.count('CAA')
            cga_count = repeat_region.count('CGA')
            cat_count = repeat_region.count('CAT')
            
            # Detect interruptions
            has_caa_interruption = caa_count > 0
            has_other_variants = (cga_count + cat_count) > 0
            
            # Calculate purity
            total_repeats = len(repeat_region) // 3
            pure_cag = (cag_count / total_repeats * 100) if total_repeats > 0 else 0
            
            # Detect fragmentation (gaps in coverage)
            # This looks at read mapping quality and cigar string
            cigar_ops = read.cigartuples
            has_deletions = any(op[0] == 2 for op in cigar_ops) if cigar_ops else False
            has_insertions = any(op[0] == 1 for op in cigar_ops) if cigar_ops else False
            is_fragmented = has_deletions or has_insertions
            
            # Classify as "clean", "interrupted", or "fragmented"
            if is_fragmented:
                structure_type = 'fragmented'
            elif has_caa_interruption or has_other_variants:
                structure_type = 'interrupted'
            else:
                structure_type = 'clean'
            
            rows.append({
                'sample': sample_name,
                'read_id': read.query_name,
                'allele': allele,
                'CAG_count': cag_count,
                'CAA_interruptions': caa_count,
                'Other_variants': cga_count + cat_count,
                'Pure_CAG_percent': pure_cag,
                'Structure_type': structure_type,
                'Has_deletions': has_deletions,
                'Has_insertions': has_insertions,
                'Read_length': read.infer_query_length(),
                'Mapping_quality': read.mapping_quality,
            })
        
        except (KeyError, TypeError):
            continue
    
    bamfile.close()
    return pd.DataFrame(rows)

# Analyze all samples
print("\n" + "="*80)
print("REPEAT STRUCTURE VARIANT ANALYSIS")
print("="*80 + "\n")

samples = {
    'Control_11': 'Control_11.trgt.sorted.spanning.bam',
    'H1': 'H1.trgt.sorted.spanning.bam',
    'H3': 'H3.trgt.sorted.spanning.bam',
}

all_structure_data = {}

for sample_name, bam_path in samples.items():
    print(f"Analyzing {sample_name}...")
    df = analyze_repeat_structure(bam_path, sample_name)
    all_structure_data[sample_name] = df
    
    # Summary statistics
    print(f"\n  Total reads: {len(df)}")
    print(f"  Clean repeats: {(df['Structure_type'] == 'clean').sum()} ({(df['Structure_type'] == 'clean').mean()*100:.1f}%)")
    print(f"  Interrupted: {(df['Structure_type'] == 'interrupted').sum()} ({(df['Structure_type'] == 'interrupted').mean()*100:.1f}%)")
    print(f"  Fragmented: {(df['Structure_type'] == 'fragmented').sum()} ({(df['Structure_type'] == 'fragmented').mean()*100:.1f}%)")
    
    # CAA interruptions
    with_caa = (df['CAA_interruptions'] > 0).sum()
    print(f"  Reads with CAA interruptions: {with_caa} ({with_caa/len(df)*100:.1f}%)")
    
    # Save
    df.to_csv(f'{sample_name}_repeat_structure.csv', index=False)
    print(f"  ✓ Saved {sample_name}_repeat_structure.csv\n")

# Create comparison of structure types
print("\n" + "="*80)
print("STRUCTURE TYPE COMPARISON")
print("="*80 + "\n")

structure_summary = []
for sample_name, df in all_structure_data.items():
    expanded = df[df['allele'] == 1]
    
    structure_summary.append({
        'Sample': sample_name,
        'Clean_%': (expanded['Structure_type'] == 'clean').mean() * 100,
        'Interrupted_%': (expanded['Structure_type'] == 'interrupted').mean() * 100,
        'Fragmented_%': (expanded['Structure_type'] == 'fragmented').mean() * 100,
        'Avg_CAA_interruptions': expanded['CAA_interruptions'].mean(),
        'Avg_purity_%': expanded['Pure_CAG_percent'].mean(),
    })

structure_df = pd.DataFrame(structure_summary)
print(structure_df.to_string(index=False))
structure_df.to_csv('HTT_structure_type_summary.csv', index=False)
print("\n✓ Saved HTT_structure_type_summary.csv")

