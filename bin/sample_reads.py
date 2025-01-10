#!/usr/bin/env python3

import argparse
import random
from Bio import SeqIO
from typing import Set
import gzip 
from functools import partial 

def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.
    
    Returns:
        argparse.Namespace: Parsed command line arguments
    """
    parser = argparse.ArgumentParser(description="Randomly sample N paired-end reads from FASTQ files")
    parser.add_argument("-f1", "--forward", required=True, help="Input forward (R1) FASTQ file")
    parser.add_argument("-f2", "--reverse", required=True, help="Input reverse (R2) FASTQ file")
    parser.add_argument("-n", "--num_reads", type=int, required=True, help="Number of reads to sample")
    parser.add_argument("-s", "--seed", type=int, default=0, help="Set the random seed to use in sampling.")
    parser.add_argument("-o1", "--out_forward", required=True, help="Output forward (R1) FASTQ file")
    parser.add_argument("-o2", "--out_reverse", required=True, help="Output reverse (R2) FASTQ file")
    return parser.parse_args()

def count_fastq_records(path: str, read_function) -> int:
    """
    Count the number of records in a FASTQ file.
    
    Args:
        filename (str): Path to FASTQ file
        
    Returns:
        int: Number of records in the file
    """
    count = 0
    with read_function(path) as f:
        for _ in SeqIO.parse(f, "fastq"):
            count += 1
    return count

def generate_random_indices(total_reads: int, num_samples: int) -> Set[int]:
    """
    Generate a set of random indices.
    
    Args:
        total_reads (int): Total number of reads in the input files
        num_samples (int): Number of reads to sample
        
    Returns:
        Set[int]: Set of random indices
    """
    if num_samples > total_reads:
        raise ValueError(f"Cannot sample {num_samples} reads from {total_reads} total reads")
    
    return set(random.sample(range(total_reads), num_samples))

def sample_fastq_files(args: argparse.Namespace) -> None:
    """
    Sample random reads from paired FASTQ files.
    
    Args:
        args (argparse.Namespace): Parsed command line arguments
    """
    random.seed(args.seed)

    read_fn = partial(gzip.open, mode='rt') if args.forward.endswith('.gz') else partial(open, mode='r')
    write_fn = partial(gzip.open, mode='wt') if args.out_forward.endswith('.gz') else partial(open, mode='w')

    # Count total reads
    total_reads = count_fastq_records(args.forward, read_fn)
    
    # Verify both files have same number of reads
    reverse_reads = count_fastq_records(args.reverse, read_fn)

    if total_reads != reverse_reads:
        raise ValueError(f"Forward and reverse files have different numbers of reads: {total_reads} vs {reverse_reads}")
    
    # Generate random indices
    random_indices = generate_random_indices(total_reads, args.num_reads)

    # Process forward reads
    with read_fn(args.forward) as in_f1, write_fn(args.out_forward) as out_f1:
        for idx, record in enumerate(SeqIO.parse(in_f1, "fastq")):
            if idx in random_indices:
                SeqIO.write(record, out_f1, "fastq")
    
    # Process reverse reads
    with read_fn(args.reverse) as in_f2, write_fn(args.out_reverse) as out_f2:
        for idx, record in enumerate(SeqIO.parse(in_f2, "fastq")):
            if idx in random_indices:
                SeqIO.write(record, out_f2, "fastq")

def main() -> None:
    """Main function to execute the script."""
    # Parse arguments
    args = parse_arguments()
    
    try:
        # Sample reads
        sample_fastq_files(args)
        print(f"Successfully sampled {args.num_reads} reads to {args.out_forward} and {args.out_reverse}")
    except Exception as e:
        print(f"Error: {str(e)}")
        exit(1)

if __name__ == "__main__":
    main()