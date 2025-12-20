import os

def get_strandedness_list():
    """Get strandedness values for all sample_runs in the same order as input files"""
    return [get_strandedness(sr) for sr in annot.index]

rule count_matrix:
    input:
        expand(
            os.path.join(result_path, "align", "{sample_run}_ReadsPerGene.out.tab"),
            sample_run=annot.index,
        ),
    output:
        os.path.join(result_path, "counts", "all_counts.tsv"),
    log:
        "logs/count_matrix.log",
    params:
        samples=annot.index.tolist(),
        strand=get_strandedness_list(),
    conda:
        "../env/pandas.yaml"
    script:
        "../scripts/count-matrix.py"
