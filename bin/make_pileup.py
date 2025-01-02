#!/usr/bin/env python3
import argparse
import pysam

def write_pileup(bamfile_path, output_path, window_size=10, min_base_quality=0, min_mapping_quality=0):
    """
    Generate pileup statistics for each position in each reference sequence from a BAM file.
    
    Args:
        bamfile_path (str): Path to the BAM file
        output_path (str): Path to output TSV file
        window_size (int): Size of sliding window
        min_base_quality (int): Minimum base quality to consider (default: 0)
        min_mapping_quality (int): Minimum mapping quality to consider (default: 0)
        plot_path (str, optional): Path to save the plot(s)
        plot_mode (str): 'facet' for single stacked plot, 'separate' for individual plots
    """
    
    # Prepare a list to collect data for DataFrame and plotting
    pileup_data = []
    
    with pysam.AlignmentFile(bamfile_path, "rb") as bamfile, open(output_path, 'w') as outfile:
        outfile.write("\t".join(["contig", "position", "depth", "A", "C", "G", "T",  "del", "N", 'base', 'max_freq', 'window']) + '\n')
        
        # Iterate through each reference sequence in the BAM file
        for reference in bamfile.references:
            sliding_window = []
            # Generate pileup for current reference sequence
            for plup_column in bamfile.pileup(
                reference,
                ignore_orphans=False,
                min_base_quality=min_base_quality,
                min_mapping_quality=min_mapping_quality
            ):
                freqs = {'A': 0, 'T': 0, 'G': 0, 'C': 0, '-': 0, 'N': 0}
                
                # Count bases at current position
                for plup_read in plup_column.pileups:
                    if plup_read.is_del:
                        freqs['-'] += 1
                    else:
                        base = plup_read.alignment.query_sequence[plup_read.query_position]
                        freqs[base] += 1
                
                # Calculate depth (excluding refskips)
                base_freqs = [freqs[base] for base in "ACGT-"]
                depth = sum(base_freqs)
                top_base = max(freqs, key=lambda x : freqs.get(x))
                max_freq = max(base_freqs) / depth if depth != 0 else 0
                max_freq = round(max_freq, 2)
            
                if len(sliding_window) > int(window_size):
                    sliding_window.pop(0)
                
                # Prepare output for TSV
                outputs = [
                    plup_column.reference_name,
                    plup_column.pos + 1,  # 1-based coordinates
                    depth,
                    freqs['A'],
                    freqs['C'],
                    freqs['G'],
                    freqs['T'],
                    freqs['-'],
                    freqs['N'], 
                    top_base,
                    max_freq,
                    "".join(sliding_window)
                ]
                outputs = [str(x) for x in outputs]
                
                outfile.write("\t".join(outputs) + '\n')
                sliding_window.append(top_base)
                
                # Collect data for DataFrame
                pileup_data.append({
                    'contig': plup_column.reference_name,
                    'position': plup_column.pos + 1,
                    'max_freq': max_freq,
                    'window': "".join(sliding_window),
                    'top_base': top_base,
                    'depth': depth
                })

 
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='BAM file of reads aligned to reference')
    parser.add_argument('-b', '--min-base-qual', type=int, default=0, help='Minimum base quality to include base in pileup. Default: 0')
    parser.add_argument('-m', '--min-map-qual', type=int, default=0, help='Minimum mapping quality to include a read in pileup. Default: 0')
    parser.add_argument('-w', '--window-size', type=int, default=10, help='Number of bases to display in the sliding window column. Default: 10')
    parser.add_argument('-o', '--output', required=True, help='Output path of pileup in TSV format')
 
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    write_pileup(
        args.input, 
        args.output, 
        args.window_size, 
        args.min_base_qual, 
        args.min_map_qual
	)
