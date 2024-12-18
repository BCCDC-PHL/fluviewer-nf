#%%
import argparse
import pandas as pd
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

DROP_COLS = "allele ANN INFO unknown transcript_biotype rank_total distance_to_feature errors_warnings_info".split()

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

def reformat_vcf(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extract snpEff annotations from INFO field and create separate columns
    
    Args:
        df (pd.DataFrame): VCF DataFrame with INFO column containing ANN field
        
    Returns:
        pd.DataFrame: DataFrame with parsed annotations
    """

    info_df = pd.DataFrame(list(df['INFO'].str.split(';').apply(lambda x: dict([y.split("=") for y in x]))))
    info_df = info_df[INFO_FIELDS]

    # Extract ANN field from INFO
    ann_df = info_df['ANN'].str.split(",").explode().str.split(r"\|", expand=True)
    ann_df.columns = SNPEFF_ANN_FIELDS

    # FOR TESTING 
    ann_df.loc[ann_df['hgvs_p'].str.contains("_"), 'hgvs_p'] = 'p.Tyr689_Gln690del'

    aa_muts = ann_df['hgvs_p'].str.extract(r'p\.([^\.]+)').squeeze()
    #aa_muts = aa_muts.loc[aa_muts.str.match(r'^[A-z]+\d+[A-z]+$',na=False)]

    mask_del_multi = aa_muts.str.match(r"[A-Za-z\d]+_[A-Za-z\d]+del$", na=False)

    aa_muts_split = aa_muts.str.extract(r'^([A-z]+)(\d+)([A-z]+)$')
    aa_muts_split = aa_muts_split.fillna("")
    
    aa_muts_split[0] = aa_muts_split[0].apply(convert_aa_three_one)
    aa_muts_split[2] = aa_muts_split[2].apply(convert_aa_three_one)


    # add logic for deletions 
    aa_muts_del = aa_muts.str.extract(r'[A-Za-z]+(\d+)_[A-Za-z]+(\d+)(del)$').fillna('')
    aa_muts_del[1] = aa_muts_del[0] + "-" + aa_muts_del[1]
    aa_muts_del[0] = ''

    aa_muts_split.loc[mask_del_multi] = aa_muts_del[mask_del_multi]
    aa_muts_split.columns = ['REF_AA', 'POS_AA', 'ALT_AA']

    ann_df_merge = pd.concat([ann_df, aa_muts_split], axis=1)
    info_df_merge = info_df.merge(ann_df_merge, how='right', left_index=True, right_index=True)
    df_final = df.merge(info_df_merge, how='right', left_index=True, right_index=True)

    def expand(ref, pos, alt):
        if pos == '':
            return ''
        ref, alt = str(ref), str(alt)
        results = []

        if alt == 'del':
            return pos + alt
        if alt == 'fs':
            return ref + pos + alt
    
        for n in range(len(ref)):
            results.append(ref[n] + str(int(pos) + n) + alt[n])

        return ";".join(results)
    
    expand_convert_aa = df_final[['REF_AA','POS_AA','ALT_AA']].apply(lambda x : expand(x['REF_AA'], x['POS_AA'], x['ALT_AA']), axis=1).to_frame()
    expand_convert_aa.columns = ['MUT_LIST_AA']
    df_final = pd.concat([df_final, expand_convert_aa], axis=1)
    
    return df_final    

def filter_vcf(df):
    df = df.drop(DROP_COLS, axis=1)
    df = df.drop(df.index[df['gene_name'].str.contains("CHR")])

    return df

def get_args():
    parser = argparse.ArgumentParser(description="Detect mutations of interest.")
    parser.add_argument("-i", "--input", required=True, help="Path to the folder containing mutation files")
    parser.add_argument("-o", "--output", required=True, help="Path to the output CSV file")

    return parser.parse_args()

def main():
    args = get_args()
    
    print("Reading VCF ...")
    vcf_df, header = read_vcf(args.input)
    
    print("Reformatting VCF ...")

    vcf_df_reformat = reformat_vcf(vcf_df)
    
    print("Filtering for relevant rows/columns ...")
    vcf_df_final = filter_vcf(vcf_df_reformat)

    
    vcf_df_final.to_csv(args.output, sep='\t', index=False)
    print("Complete.")

if __name__ == '__main__':
    main()
    


