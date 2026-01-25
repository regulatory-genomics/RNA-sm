#!/usr/bin/env python3
"""
Process RSEM gene expression files to extract:
1. Meta file: transcript_id, gene_id, length, effective_length
2. Three matrices: expected_count, TPM, FPKM (sample × gene)
"""

import sys

# logging
sys.stderr = open(snakemake.log[0], "w")

import os
import pandas as pd
import numpy as np
from pathlib import Path


def process_rsem_files(gene_files, output_files):
    """
    Process all RSEM .genes.results files and generate meta file and expression matrices.
    
    Parameters:
    -----------
    gene_files : list
        List of paths to RSEM .genes.results files
    output_files : dict
        Dictionary with keys: 'meta', 'expected_count', 'tpm', 'fpkm'
        containing output file paths
    """
    print(f"Found {len(gene_files)} gene expression files")
    
    # Initialize dictionaries to store matrices
    expected_count_dict = {}
    tpm_dict = {}
    fpkm_dict = {}
    
    # Initialize meta data list
    meta_data_list = []
    
    # Process each file
    for file_path in gene_files:
        # Extract sample name from filename
        sample_name = Path(file_path).stem.replace('.genes', '')
        print(f"Processing {sample_name}...")
        
        # Read the file
        df = pd.read_csv(file_path, sep='\t')
        
        # Extract meta information (transcript_id, gene_id, length, effective_length)
        # Note: transcript_id(s) can contain multiple IDs separated by commas
        for idx, row in df.iterrows():
            gene_id = row['gene_id']
            transcript_ids = str(row['transcript_id(s)']).split(',')
            length = row['length']
            effective_length = row['effective_length']
            
            # Expand transcript IDs - one row per transcript
            for transcript_id in transcript_ids:
                transcript_id = transcript_id.strip()
                meta_data_list.append({
                    'transcript_id': transcript_id,
                    'gene_id': gene_id,
                    'length': length,
                    'effective_length': effective_length
                })
        
        # Store expression values (using gene_id as key)
        expected_count_dict[sample_name] = df.set_index('gene_id')['expected_count']
        tpm_dict[sample_name] = df.set_index('gene_id')['TPM']
        fpkm_dict[sample_name] = df.set_index('gene_id')['FPKM']
    
    # Create meta DataFrame and remove duplicates (keep first occurrence)
    meta_df = pd.DataFrame(meta_data_list)
    meta_df = meta_df.drop_duplicates(subset=['transcript_id', 'gene_id'], keep='first')
    meta_df = meta_df.sort_values(['gene_id', 'transcript_id']).reset_index(drop=True)
    
    # Create expression matrices (sample × gene)
    # Get all unique gene IDs across all samples
    all_genes = set()
    for sample_data in [expected_count_dict, tpm_dict, fpkm_dict]:
        for sample_name, gene_data in sample_data.items():
            all_genes.update(gene_data.index)
    all_genes = sorted(list(all_genes))
    
    # Create matrices
    expected_count_matrix = pd.DataFrame(index=sorted(expected_count_dict.keys()), columns=all_genes)
    tpm_matrix = pd.DataFrame(index=sorted(tpm_dict.keys()), columns=all_genes)
    fpkm_matrix = pd.DataFrame(index=sorted(fpkm_dict.keys()), columns=all_genes)
    
    # Fill matrices
    for sample_name in expected_count_matrix.index:
        expected_count_matrix.loc[sample_name] = expected_count_dict[sample_name].reindex(all_genes, fill_value=0.0)
        tpm_matrix.loc[sample_name] = tpm_dict[sample_name].reindex(all_genes, fill_value=0.0)
        fpkm_matrix.loc[sample_name] = fpkm_dict[sample_name].reindex(all_genes, fill_value=0.0)
    
    # Save outputs
    print(f"\nSaving outputs:")
    print(f"  - Meta file: {output_files['meta']}")
    print(f"  - Expected count matrix: {output_files['expected_count']}")
    print(f"  - TPM matrix: {output_files['tpm']}")
    print(f"  - FPKM matrix: {output_files['fpkm']}")
    
    meta_df.to_csv(output_files['meta'], sep='\t', index=False)
    expected_count_matrix.to_csv(output_files['expected_count'], sep='\t')
    tpm_matrix.to_csv(output_files['tpm'], sep='\t')
    fpkm_matrix.to_csv(output_files['fpkm'], sep='\t')
    
    print(f"\nSummary:")
    print(f"  - Meta file: {len(meta_df)} transcript-gene pairs")
    print(f"  - Expression matrices: {len(all_genes)} genes × {len(expected_count_matrix)} samples")


# Get inputs and outputs from snakemake
gene_files = snakemake.input
output_files = {
    'meta': snakemake.output.meta,
    'expected_count': snakemake.output.expected_count,
    'tpm': snakemake.output.tpm,
    'fpkm': snakemake.output.fpkm
}

# Process files
process_rsem_files(gene_files, output_files)

print("\nProcessing completed successfully!")
