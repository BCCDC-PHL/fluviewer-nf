#!/usr/bin/env python3
#%%
import os 
import sys
import pandas as pd 
import argparse
from Bio import SeqIO
import traceback

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('blastn', help='BlastN results file in TSV format')
	parser.add_argument('-m', '--metric', default="bitscore", type=str, help='Scoring metric to select best reference/contig. Either bitscore (default), rawscore, or bsr.')
	parser.add_argument('-s', '--seqs', required=True, help='Sequences file in FASTA format. Either the contig file or reference BLAST database in FASTA format.')
	parser.add_argument('-o','--out_tsv', default=None, help='Filtered results output')
	parser.add_argument('-O','--out_fasta', default=None, help='Filtered FASTA output')
	parser.add_argument('--min_cov', default=25, type=int, help='Minimum coverage of database reference sequence by contig (percentage, default = 25)')
	parser.add_argument('-i','--min_id', default=90, type=int, help='Minimum nucleotide sequence identity between database reference sequence and contig (percentage, default = 90)')
	return parser.parse_args()

def load_blast_format(varname='BLAST_OUTFMT'):
	DEFAULT_OUTFMT = 'qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore score'
	blast_outfmt = os.environ.get(varname)

	if not blast_outfmt:
		print(f"BLAST_OUTFMT environment variable not found. Using default column format: {DEFAULT_OUTFMT}")
		blast_outfmt = DEFAULT_OUTFMT
	
	return blast_outfmt

def parse_blast(filepath):
	cols = load_blast_format().split()
	blast_df = pd.read_csv(filepath, sep='\t', names=cols)

	blast_df.insert(0, 'sample_name', os.path.basename(filepath).split("_")[0])

	blast_df.insert(3, 'segment', blast_df['sseqid'].str.split("|").str[0])
	return blast_df


#%%
def filter_alignments(blast_results, score_columns, min_cov, min_id):
	'''
	Find best contig for each genome segment. Returns datasheet with best contigs.'''
	print('Filtering alignments...')

	blast_results = blast_results.drop_duplicates()

	# apply a min-max scaling to the scoring columns
	scoring = blast_results[score_columns]
	scoring = (scoring - scoring.min()) / (scoring.max() - scoring.min())
	scoring['composite'] = scoring.sum(axis=1)

	blast_results = pd.concat([blast_results, scoring[['composite']]], axis=1)

	# group by contig and take only the best result
	idxmax = blast_results.groupby(['qseqid'])['composite'].idxmax()
	filtered = blast_results.loc[idxmax]

	# group by segment and take the top result
	idxmax = blast_results.groupby(['segment'])['composite'].idxmax()
	filtered = blast_results.loc[idxmax]

	filtered = filtered.sort_values('composite', ascending=False).reset_index(drop=True)

	if blast_results.shape[0] == 8:
		print("WARNING: Did not find 8 segments. Number of alignments present: ", blast_results.shape[0] )

	return filtered

def write_best_sequence(blast_results, database_path, fasta_out):
	'''
	Looks up the single top hit in the DB FASTA file and writes them to their own FASTA file.
	'''
	if blast_results.shape[0] == 0:
		print("ERROR: No BLAST results found when trying to write output sequence.", file=sys.stderr)
		sys.exit(1) 

	# Open contigs FASTA and load seqs into dict (key=seq header)
	db_seqs = SeqIO.to_dict(SeqIO.parse(database_path, 'fasta'))

	final_references = [db_seqs[ref_name] for ref_name in blast_results['sseqid']]

	# extract the top hit from the relevant sequence file 
	SeqIO.write(final_references, fasta_out, 'fasta')


def main():
	args = parse_args()

	# parse the BLAST file
	blast_df = parse_blast(args.blastn)
	
	# exit if no BLAST results are found 
	if blast_df.shape[0] == 0:
		print("WARNING: No data found in blast input file. Exiting.", file=sys.stderr)
		sys.exit(1)
	
	blast_df = filter_alignments(blast_df, args.metric.split(","), args.min_cov, args.min_id)

	# output blast results 
	blast_df.to_csv(args.out_tsv, index=False, sep='\t')

	# write the best FASTA outputs
	write_best_sequence(blast_df, args.seqs, args.out_fasta)


if __name__ == '__main__':
	main()
# %%
