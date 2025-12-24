from multiqc.utils import config

def execution_start():
    """ Code to run when MultiQC starts.
    Here we define the search pattern for our plugin files.
    """
    search_patterns = {
        'rnaseq/rnaseqqc': {
            'fn': '*metrics.tsv',  # Look for files ending in this
            # Optional: 'contents': 'Mapping Rate' (Look for specific text inside)
        },
        'rnaseq/gene_type_counts': {
            'fn': '*.json',  # Look for files ending in this (matches {sample}.json from qc_gene_type_count rule)
            'contents': 'gene_type_count'  # Ensure it's the right JSON structure
        },
        'rnaseq/rsem': {
            'fn': ['*json', "*.cnt"],
            'contents': 'num_genes_detected'
        },

    }
    config.update_dict(config.sp, search_patterns)