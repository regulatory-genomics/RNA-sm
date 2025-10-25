rule count_matrix:
    input:
        expand(
            "results/align/{unit.sample_name}-{unit.unit_name}_ReadsPerGene.out.tab",
            unit=units.itertuples(),
        ),
    output:
        "results/counts/all.tsv",
    log:
        "logs/count-matrix.log",
    params:
        samples=units["sample_name"].tolist(),
        strand=get_strandedness(units),
    script:
        "../scripts/count-matrix.py"
