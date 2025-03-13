#!/usr/bin/env python3
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
from io import StringIO
from glob import glob

script_dir = os.path.dirname(__file__)

EXPR_SEGMENT = re.compile(r"\{[^\}]+\}")
BLAST_OUTFMT6_COLS = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'.split()

def path_check(path):
	if EXPR_SEGMENT.search(path):
		return path
	else:
		raise argparse.ArgumentTypeError('ERROR: Output TSV path must include {segment} to indicate the segment variable.')

def init_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-q', '--query', required=True, help='FASTA consensus sequence containing all query sequences. Segment MUST be included in the header.')
	parser.add_argument('-d', '--database', required=True, help='FASTA reference database containing references for all segments x subtypes')
	parser.add_argument('-s', '--header_pos_seg',  default=1, help='Zero-indexed position of the segment in the consensus headers')
	parser.add_argument('-t', '--header_pos_st', default=2, help='Zero-indexed position of the subtype in the consensus headers')
	parser.add_argument('-D', '--header_delim', default='|', help='FASTA reference database containing references for all segments x subtypes')
	parser.add_argument('-o', '--output', type=path_check, help='Mutation outputs in TSV format. Requires {segment} marker to indicate the gene name in the output.') 
	parser.add_argument('-O', '--outaln', type=path_check, help='Pairwise alignments used to call mutations. Requires {segment} marker to indicate the gene name in the output.') 
	return parser

class MissingSubtypeException(Exception):
    def __init__(self, message):            
        super().__init__(message)


def make_blast_db(db_path):
	if not os.path.exists(db_path + ".pdb"):
		try: 
			command = f"makeblastdb -in {db_path} -dbtype prot".split()
			result = run(command, capture_output=True, check=True)
		except OSError as e :
			print(f"An error occurred while building the BLAST database: {str(e)}")
			raise

def run_blast(query_path, db_path):
	
	try:
		make_blast_db(db_path)
		command = f"blastx -query {query_path} -db {db_path} -outfmt 6".split()
		result = run(command, capture_output=True, check=True)

		df_raw = StringIO(result.stdout.decode('utf-8'))
		df = pd.read_csv(df_raw, sep='\t', header=None, names=BLAST_OUTFMT6_COLS)

		return df

	except Exception as e :
		print(f"An error occurred while building the BLAST database: {str(e)}")
		raise

def filter_blast_df(blast_df):
	idxmax = blast_df.groupby(['qseqid'])['bitscore'].idxmax()
	filtered = blast_df.loc[idxmax].set_index('qseqid')

	return filtered['sseqid'].to_dict()

def exit(outpath):
	mutations_df = pd.DataFrame([], columns=['CHROM','REF','POS','ALT'])
	mutations_df.to_csv(outpath, sep='\t', index=False)

#%%
def main():
	parser = init_parser()
	args = parser.parse_args()	

	# load the query sequences into a dict structure
	# sequence = qry_seqs[SEGMENT]
	qry_seqs = list(SeqIO.parse(args.query, 'fasta'))

	ref_database = SeqIO.to_dict(SeqIO.parse(args.database, 'fasta'))

	aligner = init_aligner()

	blast_df = run_blast(args.query, args.database)

	qry_align_data = filter_blast_df(blast_df)

	for nt_qry in qry_seqs:

		# something went wrong in this case; skip for now 
		if nt_qry.id not in qry_align_data:
			print("No BLAST result found for ", nt_qry.id)
			continue
		
		# retrieve the name of the most appropriate reference (based on a BLAST alignment)
		ref_name = qry_align_data[nt_qry.id]
		segment = ref_name.split("|")[0]

		# grab the appropriate reference sequence for the subtype
		aa_ref = ref_database[ref_name]

		print("Computing AA Mutations for Segment : ", segment)
		
		# aa_qry = nt_qry[nt_qry_start_position-1:].translate()  # less room for error / mistranslation with this approach
		aa_qry, _ , _ = flex_translate(nt_qry)

		# compute all mutations using the pairwise alignment
		ref_aln, qry_aln = pairwise_alignment(aligner, aa_ref, aa_qry)

		if args.outaln:
			outpath_aln = EXPR_SEGMENT.sub(segment, args.outaln)
			with open(outpath_aln, 'w') as outfile:
				outfile.write(f">{ref_name}\n{ref_aln}\n>{nt_qry.id}\n{qry_aln}\n")

		mutations = get_mutations(ref_aln, qry_aln)

		# convert list of mutations to data frame
		mutations_df = pd.DataFrame(mutations, columns=['REF','POS','ALT'])
		mutations_df.insert(0, "CHROM", aa_ref.name)

		outpath_mutations = EXPR_SEGMENT.sub(segment, args.output)

		# output mutations to CSV file 
		mutations_df.to_csv(outpath_mutations, sep='\t', index=False)

		print(f'Completed Segment: {segment}')


#%%
if __name__ == '__main__':
	main()
