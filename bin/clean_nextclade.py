#!/usr/bin/env python3
import argparse
import csv
import json
import os

def clean_tsv(input_path, output_path, config_path):
    # Read the accepted columns from the JSON file
    with open(config_path, 'r') as f:
        accepted_columns = json.load(f)['columns']

    # Read the TSV file and filter columns
    with open(input_path, 'r', newline='') as infile, open(output_path, 'w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')

        columns_present = set(accepted_columns) & set(reader.fieldnames)
        
        # Add new columns for DATASET_NAME and DATASET_VERSION
        fieldnames = accepted_columns.copy()
        fieldnames.append('DATASET_NAME')
        fieldnames.append('DATASET_VERSION')
        fieldnames.append('NEXTCLADE_VERSION')
        
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for row in reader:
            # Filter the row to only include accepted columns
            try :
                filtered_row = {}
                for col in accepted_columns:
                    if col in columns_present:
                        filtered_row[col] = row[col] 
                    else:
                        filtered_row[col] = 'N/A'

            except KeyError as e :
                print("ERROR: Issue extracting column from TSV file")
                print(e.with_traceback)
                raise 
            
            # Add DATASET_NAME and DATASET_VERSION if they are set
            filtered_row['DATASET_NAME'] = os.environ.get('NEXTCLADE_DATASET_NAME', 'UNKNOWN')
            filtered_row['DATASET_VERSION'] = os.environ.get('NEXTCLADE_DATASET_VERSION', 'UNKNOWN')
            filtered_row['NEXTCLADE_VERSION'] = os.environ.get('NEXTCLADE_VERSION', 'UNKNOWN')
            
            writer.writerow(filtered_row)

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Clean a TSV file by subsetting columns and adding dataset metadata.")
    parser.add_argument('-i', '--input', required=True, help='Path to the input TSV file.')
    parser.add_argument('-o', '--output', required=True, help='Path to save the cleaned TSV file.')
    parser.add_argument('-c', '--config', required=True, help='Path to the JSON file containing accepted columns.')

    # Parse arguments
    args = parser.parse_args()

    # Run the cleaning function
    clean_tsv(args.input, args.output, args.config)