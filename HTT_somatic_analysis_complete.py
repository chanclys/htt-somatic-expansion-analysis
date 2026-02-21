import pysam
import pandas as pd
import numpy as np

def extract_htt_reads(bam_path, sample_name):
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    rows = []
    
    for read in bamfile:
        if not read.has_tag('TR') or read.get_tag('TR') != 'HD_HTT':
            continue
        
        try:
            motif_offsets = read.get_tag('MO')
            allele_label = read.get_tag('AL')
            start_offset = read.get_tag('SO')
            end_offset = read.get_tag('EO')
            query_seq = read.query_sequence
            
            if not query_seq or not motif_offsets:
                continue
                
        except KeyError:
            continue
        
        repeat_region = query_seq[max(0, start_offset):min(len(query_seq), end_offset)]
        cag_count = len(motif_offsets)
        caa_count = repeat_region.count('CAA')
        ccg_count = repeat_region.count('CCG')
        
        classification, details = classify_repeat_structure(repeat_region, motif_offsets, cag_count, caa_count)
        
        rows.append({
            'sample': sample_name,
            'read_id': read.query_name,
            'allele': allele_label,
            'CAG_count': cag_count,
            'CAA_count': caa_count,
            'CCG_count': ccg_count,
            'repeat_length_bp': len(repeat_region),
            'classification': classification,
            'is_clean': classification == 'clean_CAG',
            'is_interrupted': classification == 'CAA_interrupted',
            'is_fragmented': classification == 'fragmented_noisy',
            'mapq': read.mapping_quality,
        })
    
    bamfile.close()
    return rows


def classify_repeat_structure(repeat_seq, motif_offsets, cag_count, caa_count):
    if caa_count > 0:
        if cag_count >= caa_count:
            return 'CAA_interrupted', f'{cag_count}CAG+{caa_count}CAA'
        else:
            return 'fragmented_noisy', 'CAA_dominant'
    
    expected_codons = len(repeat_seq) // 3
    actual_codons = cag_count
    
    if actual_codons == expected_codons:
        return 'clean_CAG', f'{cag_count}CAG'
    
    if motif_offsets:
        gaps = []
        for i in range(len(motif_offsets) - 1):
            gap = motif_offsets[i+1] - motif_offsets[i]
            if gap != 3:
                gaps.append(gap)
        
        if gaps:
            return 'frameshift_indels', f'gaps'
    
    if actual_codons < expected_codons * 0.5:
        return 'fragmented_noisy', f'only {actual_codons}/{expected_codons} codons'
    
    return 'clean_CAG', f'{cag_count}CAG'


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
    
    normal_allele_data = df[df['allele'] == 0]
    expanded_allele_data = df[df['allele'] == 1]
    
    if len(normal_allele_data) > 0:
        normal_median = normal_allele_data['CAG_count'].median()
        print(f"Normal allele median CAG: {normal_median:.0f}")
    
    if len(expanded_allele_data) > 0:
        exp_median = expanded_allele_data['CAG_count'].median()
        exp_max = expanded_allele_data['CAG_count'].max()
        print(f"Expanded allele median CAG: {exp_median:.0f}, max: {exp_max:.0f}")
    
    pct_gt_100 = (df['CAG_count'] > 100).mean() * 100
    pct_gt_500 = (df['CAG_count'] > 500).mean() * 100
    pct_gt_1000 = (df['CAG_count'] > 1000).mean() * 100
    
    print(f"% reads >100 CAG: {pct_gt_100:.2f}%")
    print(f"% reads >500 CAG: {pct_gt_500:.2f}%")
    print(f"% reads >1000 CAG: {pct_gt_1000:.2f}%")
    
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
                df.to_csv(f'{sample_name}_HTT_detailed_metrics.csv', index=False)
                print(f"Saved {sample_name}_HTT_detailed_metrics.csv\n")
        except Exception as e:
            print(f"ERROR: {e}\n")
    
    if len(all_dfs) > 1:
        comparison = []
        for sample, df in all_dfs.items():
            comparison.append({
                'Sample': sample,
                'N_reads': len(df),
                'Median_CAG': df['CAG_count'].median(),
                'Max_CAG': df['CAG_count'].max(),
                'P95_CAG': df['CAG_count'].quantile(0.95),
                'Pct_gt_100': (df['CAG_count'] > 100).mean() * 100,
                'Pct_gt_500': (df['CAG_count'] > 500).mean() * 100,
                'Pct_gt_1000': (df['CAG_count'] > 1000).mean() * 100,
            })
        comp_df = pd.DataFrame(comparison)
        comp_df.to_csv('HTT_cross_sample_comparison.csv', index=False)
        print("CROSS-SAMPLE COMPARISON:")
        print(comp_df.to_string(index=False))
        print("\nSaved HTT_cross_sample_comparison.csv")
    
    merged_df = pd.concat([df.copy() for df in all_dfs.values()], ignore_index=True)
    merged_df.to_csv('HTT_all_samples_readlevel.csv', index=False)
    print(f"\nSaved HTT_all_samples_readlevel.csv ({len(merged_df)} total reads)")


if __name__ == '__main__':
    main()
