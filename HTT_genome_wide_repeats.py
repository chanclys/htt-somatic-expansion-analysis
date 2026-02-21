import pysam
import pandas as pd
import numpy as np

def analyze_all_repeats(bam_path, sample_name):
    """Analyze ALL tandem repeats in the BAM file, not just HTT"""
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    rows = []
    
    repeat_loci = {}  # Track unique repeat loci
    
    for read in bamfile:
        if not read.has_tag('TR'):
            continue
        
        try:
            repeat_id = read.get_tag('TR')
            allele_label = read.get_tag('AL')
            start_offset = read.get_tag('SO')
            end_offset = read.get_tag('EO')
            query_seq = read.query_sequence
            
            if not query_seq:
                continue
            
            repeat_region = query_seq[max(0, start_offset):min(len(query_seq), end_offset)]
            
            if len(repeat_region) < 3:
                continue
            
            # Count repeats (generic - count most common 3bp motif)
            repeat_count = len(repeat_region) // 3
            
            if repeat_id not in repeat_loci:
                repeat_loci[repeat_id] = {'total_reads': 0, 'cag_count': 0}
            
            repeat_loci[repeat_id]['total_reads'] += 1
            repeat_loci[repeat_id]['cag_count'] += repeat_count
            
            rows.append({
                'sample': sample_name,
                'repeat_id': repeat_id,
                'read_id': read.query_name,
                'allele': allele_label,
                'repeat_length_bp': len(repeat_region),
                'repeat_count': repeat_count,
            })
        
        except KeyError:
            continue
    
    bamfile.close()
    return pd.DataFrame(rows), repeat_loci

def main():
    print("\n" + "="*60)
    print("GENOME-WIDE TANDEM REPEAT ANALYSIS")
    print("="*60)
    
    files = {
        'Control_11': 'Control_11.trgt.sorted.spanning.bam',
        'H1': 'H1.trgt.sorted.spanning.bam',
        'H3': 'H3.trgt.sorted.spanning.bam',
    }
    
    all_repeats = {}
    
    for sample_name, bam_path in files.items():
        print(f"\nAnalyzing {sample_name}...")
        df, loci = analyze_all_repeats(bam_path, sample_name)
        all_repeats[sample_name] = (df, loci)
        
        print(f"  Found {len(loci)} unique repeat loci:")
        for repeat_id, stats in sorted(loci.items(), key=lambda x: x[1]['total_reads'], reverse=True)[:10]:
            print(f"    {repeat_id}: {stats['total_reads']} reads, avg length {stats['cag_count']/stats['total_reads']:.0f} units")
    
    # Create summary table
    print("\n" + "="*60)
    print("REPEAT LOCI SUMMARY")
    print("="*60)
    
    summary_data = []
    for sample_name, (df, loci) in all_repeats.items():
        for repeat_id, stats in loci.items():
            summary_data.append({
                'Sample': sample_name,
                'Repeat_ID': repeat_id,
                'N_reads': stats['total_reads'],
                'Avg_repeat_units': stats['cag_count'] / stats['total_reads'],
            })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv('genome_wide_repeats_summary.csv', index=False)
    print(summary_df.to_string(index=False))
    print("\n✓ Saved genome_wide_repeats_summary.csv")

if __name__ == '__main__':
    main()
