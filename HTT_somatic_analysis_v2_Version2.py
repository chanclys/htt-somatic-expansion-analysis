import pysam
import pandas as pd
import numpy as np
import re

def extract_htt_reads(bam_path, sample_name):
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    rows = []
    
    for read in bamfile:
        if not read.has_tag('TR') or read.get_tag('TR') != 'HD_HTT':
            continue
        
        try:
            allele_label = read.get_tag('AL')
            start_offset = read.get_tag('SO')
            end_offset = read.get_tag('EO')
            query_seq = read.query_sequence
            
            if not query_seq:
                continue
                
        except KeyError:
            continue
        
        # Extract the repeat region from the read
        repeat_region = query_seq[max(0, start_offset):min(len(query_seq), end_offset)]
        
        if len(repeat_region) < 3:
            continue
        
        # Count CAG repeats directly from sequence
        cag_count = repeat_region.count('CAG')
        caa_count = repeat_region.count('CAA')
        ccg_count = repeat_region.count('CCG')
        
        # Classify the repeat
        classification = classify_repeat_structure(repeat_region, cag_count, caa_count)
        
        rows.append({
            'sample': sample_name,
            'read_id': read.query_name,
            'allele': allele_label,
            'CAG_count': cag_count,
            'CAA_count': caa_count,
            'CCG_count': ccg_count,
            'repeat_length_bp': len(repeat_region),
            'repeat_seq': repeat_region[:100],  # First 100bp for inspection
            'classification': classification,
            'mapq': read.mapping_quality,
        })
    
    bamfile.close()
    return rows


def classify_repeat_structure(repeat_seq, cag_count, caa_count):
    """Classify repeat structure quality"""
    
    # Calculate what percentage is CAG
    potential_codons = len(repeat_seq) // 3
    
    if potential_codons == 0:
        return 'too_short'
    
    cag_percentage = (cag_count / potential_codons) * 100 if potential_codons > 0 else 0
    
    if caa_count == 0 and cag_percentage > 95:
        return 'clean_CAG'
    elif caa_count > 0 and cag_count > caa_count:
        return 'CAA_interrupted'
    elif caa_count > 0:
        return 'CAA_dominant'
    else:
        return 'other_composition'


def analyze_sample(bam_path, sample_name):
    print(f"\n{'='*60}")
    print(f"Analyzing: {sample_name}")
    print(f"{'='*60}")
    
    rows = extract_htt_reads(bam_path, sample_name)
    df = pd.DataFrame(rows)
    
    if len(df) == 0:
        print(f"WARNING: No HTT reads found")
        return None
    
    print(f"Extracted {len(df)} HTT-spanning reads")
    
    print(f"\nCAG count statistics across all reads:")
    print(f"  Mean: {df['CAG_count'].mean():.1f}")
    print(f"  Median: {df['CAG_count'].median():.1f}")
    print(f"  Std: {df['CAG_count'].std():.1f}")
    print(f"  Min: {df['CAG_count'].min():.0f}")
    print(f"  Max: {df['CAG_count'].max():.0f}")
    print(f"  P25: {df['CAG_count'].quantile(0.25):.1f}")
    print(f"  P50: {df['CAG_count'].quantile(0.50):.1f}")
    print(f"  P75: {df['CAG_count'].quantile(0.75):.1f}")
    print(f"  P95: {df['CAG_count'].quantile(0.95):.1f}")
    print(f"  P99: {df['CAG_count'].quantile(0.99):.1f}")
    
    # By allele
    for allele in [0, 1]:
        allele_data = df[df['allele'] == allele]
        if len(allele_data) > 0:
            print(f"\nAllele {allele} (n={len(allele_data)} reads):")
            print(f"  Median CAG: {allele_data['CAG_count'].median():.1f}")
            print(f"  Max CAG: {allele_data['CAG_count'].max():.0f}")
            print(f"  P95 CAG: {allele_data['CAG_count'].quantile(0.95):.1f}")
    
    # Expansion metrics
    pct_gt_100 = (df['CAG_count'] > 100).mean() * 100
    pct_gt_500 = (df['CAG_count'] > 500).mean() * 100
    pct_gt_1000 = (df['CAG_count'] > 1000).mean() * 100
    
    print(f"\nExpansion burden:")
    print(f"  % reads >100 CAG: {pct_gt_100:.2f}%")
    print(f"  % reads >500 CAG: {pct_gt_500:.2f}%")
    print(f"  % reads >1000 CAG: {pct_gt_1000:.2f}%")
    
    # Classification
    print(f"\nRepeat structure:")
    class_counts = df['classification'].value_counts()
    for cls, count in class_counts.items():
        pct = (count / len(df)) * 100
        print(f"  {cls}: {count} ({pct:.1f}%)")
    
    return df


def main():
    files = {
        'Control_11': 'Control_11.trgt.sorted.spanning.bam',
        'H1': 'H1.trgt.sorted.spanning.bam',
        'H3': 'H3.trgt.sorted.spanning.bam',
    }
    
    all_dfs = {}
    for sample_name, bam_path in files.items():
        try:
            df = analyze_sample(bam_path, sample_name)
            if df is not None:
                all_dfs[sample_name] = df
                df.to_csv(f'{sample_name}_HTT_detailed_metrics_v2.csv', index=False)
                print(f"\nSaved {sample_name}_HTT_detailed_metrics_v2.csv")
        except Exception as e:
            print(f"ERROR: {e}")
            import traceback
            traceback.print_exc()
    
    if len(all_dfs) > 1:
        print(f"\n\n{'='*60}")
        print("CROSS-SAMPLE COMPARISON")
        print(f"{'='*60}\n")
        
        comparison = []
        for sample, df in all_dfs.items():
            comparison.append({
                'Sample': sample,
                'N_reads': len(df),
                'Mean_CAG': df['CAG_count'].mean(),
                'Median_CAG': df['CAG_count'].median(),
                'Max_CAG': df['CAG_count'].max(),
                'P95_CAG': df['CAG_count'].quantile(0.95),
                'P99_CAG': df['CAG_count'].quantile(0.99),
                'Pct_gt_100': (df['CAG_count'] > 100).mean() * 100,
                'Pct_gt_500': (df['CAG_count'] > 500).mean() * 100,
                'Pct_gt_1000': (df['CAG_count'] > 1000).mean() * 100,
            })
        comp_df = pd.DataFrame(comparison)
        comp_df.to_csv('HTT_cross_sample_comparison_v2.csv', index=False)
        print(comp_df.to_string(index=False))
        print(f"\nSaved HTT_cross_sample_comparison_v2.csv")
    
    if all_dfs:
        merged_df = pd.concat([df.copy() for df in all_dfs.values()], ignore_index=True)
        merged_df.to_csv('HTT_all_samples_readlevel_v2.csv', index=False)
        print(f"\nSaved HTT_all_samples_readlevel_v2.csv ({len(merged_df)} total reads)")


if __name__ == '__main__':
    main()