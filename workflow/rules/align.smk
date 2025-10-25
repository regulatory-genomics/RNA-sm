rule star_pe_multi:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1="results/trimmed/{sample}-{unit}_1.fq.gz",
        # paired end reads needs to be ordered so each item in the two lists match
        fq2="results/trimmed/{sample}-{unit}_2.fq.gz",  #optional
        # path to STAR reference genome index
        idx="/storage/leikaiLab/hanlitian/database/STAR_reference/star_v2",
        gtf="/storage/leikaiLab/hanlitian/database/gtf/hg38_star.gtf"
    output:
        # see STAR manual for additional output files
        aln="results/align/{sample}-{unit}_sortedByCoord.out.bam",
        log="logs/align/{sample}-{unit}_log.out",
        reads_per_gene="results/align/{sample}-{unit}_ReadsPerGene.out.tab",
        sj="results/align/{sample}-{unit}_SJ.out.tab",
        unmapped=["results/align/unmapped/{sample}-{unit}_unmapped.1.fastq.gz","results/align/unmapped/{sample}-{unit}_unmapped.2.fastq.gz"],
    log:
        "logs/align/{sample}-{unit}_star.log",
    params:
        # optional parameters
        extra=lambda wc, input: " ".join(
            [
                "--outSAMtype BAM SortedByCoordinate",
                "--quantMode GeneCounts",
                f'--sjdbGTFfile "{input.gtf}"',
            ]
        ),
    threads: 8
    wrapper:
        "v7.6.0/bio/star/align"


rule star_index:
    input:
        fasta="/storage/leikaiLab/hanlitian/database/STAR_reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa",
        gtf="/storage/leikaiLab/hanlitian/database/gtf/hg38_star.gtf",
    output:
        directory("/storage/leikaiLab/hanlitian/database/STAR_reference/star_v2"),
    log:
        "logs/star_index_genome.log",
    cache: True
    params:
        extra=lookup(within=config, dpath="params/star/index", default=""),
    threads: 4
    wrapper:
        "v7.2.0/bio/star/index"
