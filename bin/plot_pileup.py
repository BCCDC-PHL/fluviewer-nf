#!/usr/bin/env python3
import seaborn as sns
import pandas as pd 
import matplotlib.pyplot as plt
import argparse
import os 

def load_pileup(pileup_path):
    df = pd.read_csv(pileup_path, sep='\t')
    return df

def create_facet_plot(pileup_df, plot_path):
    """
    Create dot plot(s) of max_freq across genome positions for multiple contigs.
    
    Args:
        pileup_df (DataFrame): Pandas DataFrame containing pileup information
        plot_path (str): Path to save the plot(s)
        sample_name (str): Name of the BAM file for plot title
        plot_mode (str): 'facet' for single stacked plot, 'separate' for individual plots
    """
    # Convert to DataFrame
    
    # Get unique contigs
    contigs = pileup_df['contig'].unique()
    
    plt.figure(figsize=(20, 4 * len(contigs)))  # Adjust height based on number of contigs
        
    for i, contig in enumerate(contigs, 1):
        plt.subplot(len(contigs), 1, i)
        
        contig_data = pileup_df[pileup_df['contig'] == contig]
        
        # Separate normal and low-frequency points
        normal_points = contig_data[contig_data['max_freq'] >= 0.9]
        low_freq_points = contig_data[contig_data['max_freq'] < 0.9]
        
        # Plot normal points
        plt.scatter(normal_points['position'], normal_points['max_freq'], 
                    alpha=0.3, color='gray', s=10)
        
        # Plot and label low-frequency points with annotations
        for _, row in low_freq_points.iterrows():
            plt.scatter(row['position'], row['max_freq'], 
                        color='red', s=50, edgecolors='black', linewidth=1, zorder=5)
            plt.annotate(
                f"Pos: {row['position']}\nFreq: {row['max_freq']:.2f}\nDepth: {row['depth']}\nWindow: {row['window']}",
                (row['position'], row['max_freq']),
                xytext=(10, 10),
                textcoords='offset points',
                fontsize=8,
                bbox=dict(boxstyle="round,pad=0.3", fc="yellow", ec="black", lw=1, alpha=0.7),
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0.2")
            )
        
        plt.title(f'Max Base Frequency - {contig}', fontsize=10)
        plt.xlabel('Genome Position')
        plt.ylabel('Max Base Frequency')
        plt.ylim(0.2, 1.02)
        plt.grid(True, linestyle='--', alpha=0.7)
        
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()

def create_separate_plots(pileup_df, plot_path, sample_name):

    def extract_contig_name(contig):
        """
        Extract contig name between first and second '|' characters.
        If '|' is not found twice, return the original contig name.
        """
        parts = contig.split('|')
        return parts[1] if len(parts) >= 2 else contig

    # Create separate plots for each contig
    contigs = pileup_df['contig'].unique()

    for contig in contigs:
        print(contig)
        plt.figure(figsize=(20, 6))
        
        contig_data = pileup_df[pileup_df['contig'] == contig]
        
        # Separate normal and low-frequency points
        normal_points = contig_data[contig_data['max_freq'] >= 0.9]
        low_freq_points = contig_data[contig_data['max_freq'] < 0.9]
        
        # Plot normal points
        plt.scatter(normal_points['position'], normal_points['max_freq'], 
                    alpha=0.3, color='gray', s=10)
        
        # Plot and label low-frequency points with annotations
        for _, row in low_freq_points.iterrows():
            plt.scatter(row['position'], row['max_freq'], 
                        color='red', s=50, edgecolors='black', linewidth=1, zorder=5)
            plt.annotate(
                f"Pos: {row['position']}\nFreq: {row['max_freq']:.2f}\nDepth: {row['depth']}\nWindow: {row['window']}",
                (row['position'], row['max_freq']),
                xytext=(10, 10),
                textcoords='offset points',
                fontsize=8,
                bbox=dict(boxstyle="round,pad=0.3", fc="yellow", ec="black", lw=1, alpha=0.7),
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0.2")
            )
        
        # Extract simplified contig name
        simple_contig_name = extract_contig_name(contig)
        
        plt.title(f'Max Base Frequency Across {simple_contig_name}\n{sample_name}', fontsize=16)
        plt.xlabel('Genome Position', fontsize=12)
        plt.ylabel('Max Base Frequency', fontsize=12)
        plt.ylim(0.2, 1.02)
        plt.grid(True, linestyle='--', alpha=0.7)
        
        # Save individual contig plot with unique filename using simplified contig name
        individual_plot_path = plot_path.replace('.png', f'_{simple_contig_name}.png')
        plt.tight_layout()
        plt.savefig(individual_plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='Pileup in TSV format')
    parser.add_argument('-o', '--output', required=True, help='Output path for dot plot PNG')
    parser.add_argument('-m', '--plot-mode', choices=['facet', 'separate'], default='facet', 
                        help='Plot mode: "facet" for stacked plot, "separate" for individual contig plots. Default: facet')
 
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    sample_name = os.path.basename(args.input).split("_")[0]

    pileup_df = load_pileup(args.input)

    if args.plot_mode == 'facet':
        print("Generating facet plot...")
        create_facet_plot(pileup_df, args.output)
        
    else:
        print("Generating separate plots...")
        create_separate_plots(pileup_df, args.output, sample_name)


