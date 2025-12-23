import os

rule count_matrix:
    input:
        expand(
            os.path.join(result_path, "downstream_res", "counts", "{sample}_ReadsPerGene.out.tab"),
            sample=samples.keys(),
        ),
    output:
        os.path.join(result_path, "downstream_res", "counts", "all_counts.tsv"),
    log:
        os.path.join(result_path, "logs", "count_matrix.log"),
    params:
        samples=list(samples.keys()),
        strand=get_strandedness_list(),
    conda:
        "../env/pandas.yaml"
    script:
        "../scripts/count-matrix.py"
