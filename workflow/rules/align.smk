import os

rule star_align:
    input:
        fq1=lambda w: os.path.join(result_path, "middle_file", "trimmed", f"{w.sample_run}_1.fq.gz"),
        fq2=lambda w: os.path.join(result_path, "middle_file", "trimmed", f"{w.sample_run}_2.fq.gz"),
        idx=config["resources"]["star_index"],
        gtf=config["resources"]["gtf"]
    output:
        aln=os.path.join(result_path, "Important_processed", "Bam", "{sample_run}_sortedByCoord.out.bam"),
        log_final=os.path.join(result_path, "Report", "star", "{sample_run}_Log.final.out"),
        reads_per_gene=os.path.join(result_path, "downstream_res", "counts", "{sample_run}_ReadsPerGene.out.tab"),
        sj=os.path.join(result_path, "Important_processed", "Bam", "{sample_run}_SJ.out.tab"),
        unmapped=[
            os.path.join(result_path, "middle_file", "unmapped", "{sample_run}_unmapped.1.fastq.gz"),
            os.path.join(result_path, "middle_file", "unmapped", "{sample_run}_unmapped.2.fastq.gz")
        ],
    log:
        os.path.join(result_path, "logs", "star", "{sample_run}.log"),
    params:
        extra=lambda wc, input: " ".join([
            "--outSAMtype BAM SortedByCoordinate",
            "--quantMode GeneCounts",
            f'--sjdbGTFfile "{input.gtf}"',
        ]),
        index=lambda wc, input: input.idx,
    threads: 8
    resources:
        mem_mb=50000,
        runtime=400,
    wrapper:
        "v7.6.0/bio/star/align"

rule merge_sample_runs:
    input:
        lambda w: expand(
            os.path.join(result_path, "Important_processed", "Bam", "{sample_run}_sortedByCoord.out.bam"),
            sample_run=get_runs_for_sample(w.sample)
        )
    output:
        os.path.join(result_path, "Important_processed", "Bam", "{sample}.bam")
    log:
        os.path.join(result_path, "logs", "samtools", "merge", "{sample}.log")
    threads: 4
    wrapper:
        "v3.3.3/bio/samtools/merge"

rule samtools_index_bam:
    input:
        os.path.join(result_path, "Important_processed", "Bam", "{sample_run}_sortedByCoord.out.bam")
    output:
        os.path.join(result_path, "Important_processed", "Bam", "{sample_run}_sortedByCoord.out.bam.bai")
    log:
        os.path.join(result_path, "logs", "samtools", "index", "{sample_run}.log")
    threads: 2
    wrapper:
        "v3.3.3/bio/samtools/index"

rule samtools_index_merged:
    input:
        os.path.join(result_path, "Important_processed", "Bam", "{sample}.bam")
    output:
        os.path.join(result_path, "Important_processed", "Bam", "{sample}.bam.bai")
    log:
        os.path.join(result_path, "logs", "samtools", "index", "{sample}_merged.log")
    threads: 2
    wrapper:
        "v3.3.3/bio/samtools/index"

rule star_index:
    input:
        fasta=config["resources"]["fasta"],
        gtf=config["resources"]["gtf"],
    output:
        directory(config["resources"]["star_index"]),
    log:
        os.path.join(result_path, "logs", "star_index.log"),
    cache: True
    params:
        extra=config["params"]["star"]["index"]["extra"],
    threads: 8
    resources:
        mem_mb=50000,
        runtime=100,
    wrapper:
        "v7.2.0/bio/star/index"
