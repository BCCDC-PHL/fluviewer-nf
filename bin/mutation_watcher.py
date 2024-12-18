#!/usr/bin/env python3
#%%
import argparse
import pandas as pd
from glob import glob
import os
import re 
from typing import List

#%%
SEGMENTS = 'PB2 PB1 PA HA NP NA M NS'.split()
SEGMENT_REGEX = re.compile(r'_([^_]+)_mutations')

def load_watchlist(watchlist_path):
    """
    Load the watchlist CSV file using pandas and convert it to a dict.
    """
    try:
        watchlist_df = pd.read_csv(watchlist_path)
        # Assuming the mutations are in a column named 'mutation'
        # Modify this according to your CSV structure
        watchlist_dict = watchlist_df.groupby('segment')['mutation'].apply(set).to_dict()
        return watchlist_dict
    except Exception as e:
        print(f"Error loading watchlist: {str(e)}")
        return dict()
    
def load_vcf(vcf_path):
    df = pd.read_csv(vcf_path, sep='\t', keep_default_na=False)
    mutation_dict = df.groupby('gene_name')['MUT_LIST_AA'].apply(lambda x : ";".join([y for y in x if y])).str.split(";").to_dict()
    mutation_dict = {gene : set(muts) for gene, muts in mutation_dict.items()}
    return mutation_dict

def detect_mutations(mutation_dict, watchlist_dict):
    """
    Scan the list of pandas dataframes to detect mutations in the watchlist.
    """

    detected_mutations = {}
    for segment, watch_muts in watchlist_dict.items():

        if segment in mutation_dict:
            flagged = mutation_dict[segment].intersection(watch_muts)
            if flagged:
                flagged = sorted(flagged, key = lambda x : int(re.search(r'(\d+)', x).group(1)))
                detected_mutations[segment] = flagged
    return detected_mutations

def parse_args():
    parser = argparse.ArgumentParser(description="Detect mutations of interest.")
    parser.add_argument("-i", "--input", required=True, help="Path to the folder containing mutation files")
    parser.add_argument("-w", "--watchlist", required=True, help="Path to the CSV file containing the watchlist")
    parser.add_argument("-o", "--output", required=True, help="Path to the output CSV file")
    parser.add_argument("--glob_pattern", default="*mutations.tsv", help="Glob pattern for input files (default: *.csv)")

    return parser.parse_args()
#%%
def main():

    args = parse_args()

    # Load data files
    mutation_dict = load_vcf(args.input)
    sample_name = args.input.split("_")[0]
    print(f"Loaded VCF.")

    # Load watchlist
    watchlist = load_watchlist(args.watchlist)

    # Detect mutations
    detected_mutations = detect_mutations(mutation_dict, watchlist)

    # Print results
    if detected_mutations:
        print("Detected mutations of interest:")
        dfs = []
        for segment, mutations in detected_mutations.items():
            dfs.append(pd.DataFrame({'name': sample_name, 'segment': segment, 'mutations':mutations}))
        final_df = pd.concat(dfs)
    else:
        print("No mutations of interest detected.")
        final_df = pd.DataFrame({"name":[], 'segment':[], 'mutations':[]})

    final_df.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()
