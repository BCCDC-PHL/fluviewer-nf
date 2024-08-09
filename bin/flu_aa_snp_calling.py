#!/usr/bin/env python3
"""
This flu specific script performs the following:
- translates nucleotide sequence queries into amino acid format
- searches nucleotide sequences against an amino acid sequence database using BLASTX to find best reference
- computes amino acid mutations relative to best chosen reference 
"""


#%%
import argparse
import os, sys, re
import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from tools import flex_translate, init_aligner, pairwise_alignment, get_mutations
import yaml
from subprocess import run, PIPE
from glob import glob
import traceback 

script_dir = os.path.dirname(__file__)

EXPR_SEGMENT = re.compile(r"\{[^\}]+\}")

def path_check(path):
    if EXPR_SEGMENT.search(path):
        return path
    else:
        raise argparse.ArgumentTypeError('ERROR: Output TSV path must include {segment} to indicate the segment variable.')

def init_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--query', required=True, help='Nucleotide sequences in FASTA format. Segment MUST be included in the header.')
    parser.add_argument('-d', '--db', required=True, help='Amino acid reference database in FASTA format containing all 8 segment sequences')
    parser.add_argument('-o', '--output', type=path_check, help='Mutation outputs in TSV format. Requires {segment} marker to indicate the gene name in the output.') 
    parser.add_argument('-O', '--outseq', help='Output path of translated amino acid sequences.') 
    return parser

class MissingSubtypeException(Exception):
    def __init__(self, message):            
        super().__init__(message)


def run_blastx(consensus_path, ref_db_path):

    sample_name = next(SeqIO.parse(consensus_path, 'fasta')).id.split("|")[0]
    blast_out_path = f"{sample_name}_blastx.tsv"

    try: 
        if not os.path.exists(ref_db_path + ".pdb"):
            command_args = [
                'makeblastdb', '-in', ref_db_path, '-dbtype' , 'nucl' 
            ]

            result = run(command_args, check_output=True)

    except Exception as e:
        print("ERROR: Failed to build BLASTX database")
        print(traceback.format_exc())
        return None

    try:
        command_args = [
            'blastx', '-query', consensus_path, '-db', ref_db_path, '-outfmt', '6'
        ]
        with open(blast_out_path, 'w') as outfile:
            result = run(command_args, stdout=outfile, stderr=PIPE)

    except Exception as e:
        print("ERROR: Failed to run BLASTX search.")
        print(traceback.format_exc())
        return None

    return blast_out_path


def make_reference_dict(blast_path):
    blast6cols = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

    blastdf = pd.read_csv(blast_path, sep='\t', names=blast6cols.split())
    idxmax = blastdf.groupby(['qseqid'])['bitscore'].idxmax()
    filtered = blastdf.loc[idxmax]

    return filtered.set_index("qseqid")['sseqid'].to_dict()

#%%
def main():
    parser = init_parser()
    args = parser.parse_args()	

    # load in the query consensus sequences 
    qry_seqs = list(SeqIO.parse(args.query, 'fasta'))

    # load in the reference database sequences 
    ref_database = SeqIO.to_dict(SeqIO.parse(args.db, 'fasta'))

    # run blastx on the query sequences 
    blast_result = run_blastx(args.query, args.db)

    if not blast_result:
        print("ERROR: run_blastx() failed")
        exit(1)

    ref_dict = make_reference_dict(blast_result)

    aligner = init_aligner()

    aa_qry_seqs = [] 

    for qry_seq in qry_seqs:

        # something went wrong in this case; skip for now 
        if qry_seq.id not in ref_dict:
            print("No BLAST result found for ", qry_seq.id)
            continue
        
        # retrieve the name of the most appropriate reference (based on a BLAST alignment)
        ref_name = ref_dict[qry_seq.id]
        
        # extract the segment from the header
        segment = ref_name.split("|")[0]

        # translate the nucleotide query into amino acid 
        aa_qry, _ , _ = flex_translate(qry_seq)

        # save the aa_qry sequences for output 
        aa_qry_seqs.append(aa_qry)

        # grab the appropriate reference sequence for the subtype
        aa_ref = ref_database[ref_name]

        # compute all mutations using the pairwise alignment
        ref_aln, qry_aln = pairwise_alignment(aligner, aa_ref, aa_qry)

        mutations = get_mutations(ref_aln, qry_aln)

        # convert list of mutations to data frame
        mutations_df = pd.DataFrame(mutations, columns=['POS','REF','ALT'])
        mutations_df.insert(0, "CHROM", aa_ref.name)

        outpath = EXPR_SEGMENT.sub(segment, args.output)

        # output mutations to CSV file t
        mutations_df.to_csv(outpath, sep='\t', index=False)

        print(f'Completed AA SNP Calling for Segment {segment}')

    if args.outseq:
        SeqIO.write(aa_qry_seqs, args.outseq, 'fasta')


if __name__ == '__main__':
    main()
# %%
