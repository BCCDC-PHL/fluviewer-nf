#!/usr/bin/env python3
#%%
import argparse
import pandas as pd, numpy as np
from glob import glob
import os
import re 
from tools import read_vcf
#%%
SEGMENTS = 'PB2 PB1 PA HA NP NA M NS'.split()
SEGMENT_REGEX = re.compile(r'_([^_]+)_mutations')

SNPEFF_ANN_FIELDS = [
    "allele",
    "annotation",
    "putative_impact",
    "gene_name",
    "gene_id",
    "feature_type",
    "feature_id",
    "transcript_biotype",
    "rank_total",
    "hgvs_c",
    "hgvs_p",
    "cdna_position_cdna_len",
    "cds_position_cds_len",
    "protein_position_protein_len",
    "distance_to_feature",
    "errors_warnings_info"
]

AA_CODES = {
    'Ala': 'A',
    'Arg': 'R',
    'Asn': 'N',
    'Asp': 'D',
    'Asx': 'B',
    'Cys': 'C',
    'Glu': 'E',
    'Gln': 'Q',
    'Glx': 'Z',
    'Gly': 'G',
    'His': 'H',
    'Ile': 'I',
    'Leu': 'L',
    'Lys': 'K',
    'Met': 'M',
    'Phe': 'F',
    'Pro': 'P',
    'Ser': 'S',
    'Thr': 'T',
    'Trp': 'W',
    'Tyr': 'Y',
    'Val': 'V',
    'del': 'del',
    'Ter': '*'
}

INFO_FIELDS = "DP RO AO AF CIGAR TYPE ANN".split()
DROP_COLS = "ID allele FILTER ANN INFO unknown transcript_biotype rank_total distance_to_feature errors_warnings_info".split()

def convert_aa_three_one(aa_three_str):
    aa_three_str = str(aa_three_str)
    if aa_three_str in ['', 'nan']:
        return ''
    if aa_three_str == 'fs':
        return 'fs'
    if len(aa_three_str) % 3 != 0:
        print(f'ERROR: AA notation not multiple of 3: {aa_three_str}')
        return None

    aa_one_str = ''
    for i in range(0, len(aa_three_str), 3):
        aa = aa_three_str[i:i+3]
        if aa not in AA_CODES:
            print(f'ERROR: AA notation {aa} not found within string {aa_three_str}')
            return None

        aa_one_str += AA_CODES[aa]
    return aa_one_str

def extract_info_df(df):
    # extract INFO fields from the VCF 
    info_df = pd.DataFrame(list(df['INFO'].str.split(';').apply(lambda x: dict([y.split("=") for y in x]))))
    info_df = info_df[INFO_FIELDS]
    #info_df['DP'] = info_df['DP'].astype(int)
    #info_df['VAF'] = info_df.apply(lambda row : ",".join([int(x) / row['DP'] if row['DP'] != 0 else 0 for x in row['AO'].split(",")]) , axis=1)
    return info_df

def extract_ann_df(info_df):
    # Extract ANN fields from INFO dataframe
    ann_df = info_df['ANN'].str.split(",").explode().str.split(r"\|", expand=True)
 
    # rename columns
    ann_df.columns = SNPEFF_ANN_FIELDS

    aa_muts = ann_df['hgvs_p'].str.extract(r'p\.([^\.]+)').squeeze()
    #aa_muts = aa_muts.loc[aa_muts.str.match(r'^[A-z]+\d+[A-z]+$',na=False)]

    # mask to search for multi-base deletions 
    mask_del_multi = aa_muts.str.match(r"[A-Za-z\d]+_[A-Za-z\d]+del$", na=False)

    # splits the amino acid mutations into 3 columns: ref, pos, alt 
    aa_muts_split = aa_muts.str.extract(r'^([A-z]+)(\d+)([A-z]+)$')
    aa_muts_split = aa_muts_split.fillna("")
    
    # apply the conversion function to change ref and alt into single letter amino acids 
    aa_muts_split[0] = aa_muts_split[0].apply(convert_aa_three_one)
    aa_muts_split[2] = aa_muts_split[2].apply(convert_aa_three_one)

    # extract deletion entries and reformat them into an appropriate 
    aa_muts_del = aa_muts.str.extract(r'[A-Za-z]+(\d+)_[A-Za-z]+(\d+)(del)$').fillna('')
    aa_muts_del[1] = aa_muts_del[0] + "-" + aa_muts_del[1]
    aa_muts_del[0] = ''

    aa_muts_split.loc[mask_del_multi] = aa_muts_del[mask_del_multi]
    aa_muts_split.columns = ['REF_AA', 'POS_AA', 'ALT_AA']

    return pd.concat([ann_df, aa_muts_split], axis=1)

def split_aa_muts(df):
    def split_mutations(ref, pos, alt):
        if pos == '':
            return np.nan
        ref, alt = str(ref), str(alt)
        results = []

        if alt == 'del':
            return pos + alt
        if alt == 'fs':
            return ref + pos + alt
    
        for n in range(len(ref)):
            results.append(ref[n] + str(int(pos) + n) + alt[n])

        return ";".join(results)
    
    split_aa_df = df[['REF_AA','POS_AA','ALT_AA']].apply(lambda x : split_mutations(x['REF_AA'], x['POS_AA'], x['ALT_AA']), axis=1).to_frame()
    split_aa_df.columns = ['MUT_LIST_AA']

    return split_aa_df


def reformat_vcf(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extract snpEff annotations from INFO field and create separate columns
    
    Args:
        df (pd.DataFrame): VCF DataFrame with INFO column containing ANN field
        
    Returns:
        pd.DataFrame: DataFrame with parsed annotations
    """
    # extract info and annotation dataframes 
    info_df = extract_info_df(df)
    ann_df = extract_ann_df(info_df)

    # merge dataframes back into original
    info_df_merge = info_df.merge(ann_df, how='right', left_index=True, right_index=True)
    df_merge = df.merge(info_df_merge, how='right', left_index=True, right_index=True)

    # split up the amino acid entries 
    split_aa_df = split_aa_muts(df_merge)

    # build the final nucleotide dataframe
    df_nt_final = pd.concat([df_merge, split_aa_df], axis=1)

    # # reset indices 
    # df_merge = df_merge.reset_index(drop=True)
    # split_aa_df = split_aa_df.reset_index(drop=True)

    # # explode the amino acid list values to create an AA-focused dataframe 
    # split_aa_df['MUT_LIST_AA'] = split_aa_df['MUT_LIST_AA'].str.split(';')
    # split_aa_df_explode = split_aa_df.explode('MUT_LIST_AA').reset_index(names='MUT_INDEX')

    # # build final AA-focused dataframe
    # df_aa_final = df_merge.merge(split_aa_df_explode, how='right', left_index=True, right_on='MUT_INDEX').set_index('MUT_INDEX')

    # reorganize columns
    def shift_column(df: pd.DataFrame, col_name: str, pos: int):
        df.insert(pos, col_name, df.pop(col_name))

    for col_name in reversed(['REF_AA', 'POS_AA', 'ALT_AA', 'MUT_LIST_AA']):
        shift_column(df_nt_final, col_name, 6)
        #shift_column(df_aa_final, col_name, 6)

    return df_nt_final #, df_aa_final 

def filter_vcf(df):
    df = df.drop(DROP_COLS, axis=1)
    df = df.drop(df.index[df['gene_name'].str.contains("CHR")])
    df = df.drop(df.index[df['hgvs_p']==''])
    return df

def validate_outpath(outpath):
    if not outpath.endswith('.tsv') and not outpath.endswith('.csv'):
        raise ValueError("ERROR: Output path must end with either .tsv or .csv")
    return outpath

def get_args():
    parser = argparse.ArgumentParser(description="Detect mutations of interest.")
    parser.add_argument("-i", "--input", required=True, help="Path to the folder containing mutation files")
    parser.add_argument("-o", "--output", required=True, type=validate_outpath, help="Path to the output CSV file")

    return parser.parse_args()

def main():
    args = get_args()
    
    print("Reading VCF ...")
    vcf_df, _ = read_vcf(args.input)
    
    print("Reformatting VCF ...")

    vcf_df = reformat_vcf(vcf_df)
    
    print("Filtering for relevant rows/columns ...")
    vcf_df_final = filter_vcf(vcf_df)
    
    vcf_df_final.to_csv(args.output, sep='\t', index=False)
    print("Complete.")

if __name__ == '__main__':
    main()
    



# %%
