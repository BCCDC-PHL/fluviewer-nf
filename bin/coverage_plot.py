#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict
import logging
import argparse

def get_args() -> argparse.Namespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Create coverage plots from samtools depth output"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input TSV file from samtools depth -a",
        type=Path
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output directory for coverage plots",
        type=Path
    )
    return parser.parse_args()

def read_coverage_data(input_file: Path) -> Dict[str, pd.DataFrame]:
    """Read the TSV file and organize data by segment"""
    try:
        logging.info(f"Reading coverage data from {input_file}")
        
        # Read TSV file with specified column names
        df = pd.read_csv(
            input_file,
            sep='\t',
            names=['segment', 'position', 'depth'],
            dtype={
                'segment': str,
                'position': int,
                'depth': int
            }
        )
        
        # Store grouped DataFrames directly
        coverage_data = {seg_name.split("|")[0] : sub_df for seg_name, sub_df in df.groupby('segment')}
        
        logging.info(f"Successfully processed {len(coverage_data)} segments")
        return coverage_data
        
    except Exception as e:
        logging.error(f"Error reading coverage data: {str(e)}")
        raise

def create_coverage_plots(coverage_data: Dict[str, pd.DataFrame], output_path: Path) -> None:
    """Create a single figure with subplots for all segments"""
    try:
        # Define desired segment order
        segment_order = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'M', 'NS']
        
        # Create figure with subplots
        fig, axes = plt.subplots(nrows=8, ncols=1, figsize=(12, 24))
        
        # Remove the extra space at the top by adjusting layout
        plt.subplots_adjust(hspace=0.4, top=0.98)
        
        # Create plots for each segment in specified order
        for n, segment in enumerate(segment_order):

            if segment in coverage_data:
                segment_df = coverage_data[segment]
                ax = axes[n]
                
                # Plot coverage with log scale
                ax.plot(segment_df['position'], segment_df['depth'])
                ax.set_yscale('log')
                ax.set_ylim(1, 5000)
                ax.set_title(f"Segment: {segment}")
                ax.set_xlabel("Position")
                ax.set_ylabel("Depth (log10)")
                ax.grid(True)
            else:
                # If segm   ent not found, create empty plot
                ax = axes[n]
                ax.set_title(f"Segment: {segment} (Not Found)")
                ax.set_xlabel("Position")
                ax.set_ylim(1, 5000)
                ax.set_ylabel("Depth (log10)")
                ax.grid(True)
        
        # Save plot
        
        plt.savefig(output_path, bbox_inches='tight', dpi=300)
        plt.close()
        
        logging.info(f"Created combined coverage plot at {output_path}")
        
    except Exception as e:
        logging.error(f"Error creating coverage plots: {str(e)}")
        raise

def main():
    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # Get command line arguments
    args = get_args()
    
    ## Process and create plots for all segments
    try:
        # Read data
        coverage_data = read_coverage_data(args.input)
        
        # Create combined plot
        create_coverage_plots(coverage_data, args.output)
            
        logging.info("Completed processing all segments")
        
    except Exception as e:
        logging.error(f"Error in processing: {str(e)}")
        raise

if __name__ == "__main__":
    main()
