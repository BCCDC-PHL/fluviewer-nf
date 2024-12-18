#%%
import os, sys, re
import numpy as np
import pandas as pd 
from typing import Tuple, List
from functools import partial
import gzip 

BASE_DICT = dict(zip(list("ACGT"), range(4)))


SEGMENTS = 'PB2 PB1 PA HA NP NA M NS'.split()

def read_vcf(file_path: str , max_header_lines=300) -> Tuple[pd.DataFrame, List[str]]:
    """
    Parse a Variant Call Format (VCF) file into a pandas DataFrame and header list.
    
    Args:
        file_path: Path to the VCF file
        
    Returns:
        Tuple containing:
            - pandas DataFrame with VCF data
            - List of header lines (including metadata)
    """
    # Store header lines
    header_lines = []
    
    # Read header lines
    read_function = partial(open, mode='r') if not file_path.endswith('.gz') else partial(gzip.open, mode='rt')

    with read_function(file_path) as infile:
        line_count = 0 
        line = infile.readline()
        while not line.lstrip('#').startswith("CHROM"):
            line_count += 1
            header_lines.append(line.strip())
            line = infile.readline()

            if line_count > max_header_lines:
                print(f"ERROR: Exceeded {max_header_lines} lines to search for the CHROM start point.")
                sys.exit(1)

        header_lines.append(line.strip())
        
        column_names = header_lines[-1].lstrip("#").split("\t")

        df = pd.read_csv(infile, sep='\t', names=column_names)
    
    return df, header_lines


def write_vcf(df, header_lines, output_path):
    """
    Write out a Variant Call Format (VCF) file using a Pandas df, list of headers, and output path
    
    Args:
        file_path: Path to the VCF file
        
    Returns:
        Tuple containing:
            - pandas DataFrame with VCF data
            - List of header lines (including metadata)
    """

    # Read header lines
    with open(output_path, 'w') as outfile:
        for line in header_lines:
            outfile.write(line + '\n')

        df.to_csv(outfile, sep='\t', index=False, header=False)