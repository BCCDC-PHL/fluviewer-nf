#!/usr/bin/env python3
import argparse
import pysam
import sys

def write_pileup(bamfile_path, output_path, min_base_quality=0, min_mapping_quality=0):
    """
    Generate pileup statistics for each position in each reference sequence from a BAM file.
    
    Args:
        bamfile_path (str): Path to the BAM file
        min_base_quality (int): Minimum base quality to consider (default: 0)
        min_mapping_quality (int): Minimum mapping quality to consider (default: 0)
    """
    
    # Print header
    
    with pysam.AlignmentFile(bamfile_path, "rb") as bamfile, open(output_path, 'w') as outfile:

        outfile.write("\t".join(["contig", "position", "depth", "A", "C", "G", "T",  "del", "N", 'base', 'window']) + '\n')

        # Iterate through each reference sequence in the BAM file
        for reference in bamfile.references:
            # Generate pileup for current reference sequence

            sliding_window = []

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
                max_freq = round(max_freq,2)

                sliding_window.append(top_base)

                if len(sliding_window) > 10:
                    sliding_window.pop(0)

                # Prepare output
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
                    "".join(sliding_window)
                ]

                outputs = [str(x) for x in outputs]
                
                outfile.write("\t".join(outputs) + '\n')
            


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='BAM file of reads aligned to reference' )
    parser.add_argument('-b', '--min-base-qual', default=0, help='Minimum base quality to include base in pileup. Default: 0' )
    parser.add_argument('-m', '--min-map-qual', default=0, help='Minimum mapping quality to include a read in pileup. Default: 0' )
    parser.add_argument('-o', '--output', required=True, help='Output path of pileup in TSV format' )
 
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    write_pileup(args.input, args.output, args.min_base_qual, args.min_map_qual)
