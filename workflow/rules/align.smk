import os

def get_all_trimmed_fastqs(sample_name):
    """
    Get all trimmed FASTQ files for all runs of a sample.
    Returns (r1_list, r2_list) where each list contains paths for all runs.
    """
    sample_runs = get_runs_for_sample(sample_name)
    r1_files = []
    r2_files = []
    
    for sample_run in sample_runs:
        r1_path = os.path.join(result_path, "middle_file", "trimmed", f"{sample_run}_1.fq.gz")
        r2_path = os.path.join(result_path, "middle_file", "trimmed", f"{sample_run}_2.fq.gz")
        r1_files.append(r1_path)
        r2_files.append(r2_path)
    
    return r1_files, r2_files

rule star_align:
    """
    Align all runs for a sample together using STAR.
    Takes all trimmed FASTQ files from all runs of a sample as input.
    """
    input:
        fq1=lambda w: get_all_trimmed_fastqs(w.sample)[0],
        fq2=lambda w: get_all_trimmed_fastqs(w.sample)[1],
        idx=config["resources"]["star_index"],
        gtf=config["resources"]["gtf"]
    output:
        aln=os.path.join(result_path, "Important_processed", "Bam", "{sample}_sortedByCoord.out.bam"),
        transcriptome_bam=os.path.join(result_path, "Important_processed", "Bam", "{sample}_Aligned.toTranscriptome.out.bam"),
        reads_per_gene=os.path.join(result_path, "downstream_res", "counts", "{sample}_ReadsPerGene.out.tab"),
        sj=os.path.join(result_path, "Important_processed", "Bam", "{sample}_SJ.out.tab"),
        log_final=os.path.join(result_path, "Report", "star", "{sample}_Log.final.out"),
        # STAR outputs unmapped as plain text, we will gzip them manually
        unmapped=[
            os.path.join(result_path, "middle_file", "unmapped", "{sample}_unmapped.1.fastq.gz"),
            os.path.join(result_path, "middle_file", "unmapped", "{sample}_unmapped.2.fastq.gz")
        ]
    log:
        os.path.join(result_path, "logs", "star", "{sample}.log")
    params:
        # Define a temporary prefix for this sample to avoid collisions
        prefix=lambda w: os.path.join(result_path, "middle_file", "temp_star", w.sample, ""),
        # Join all R1 and R2 files with commas for STAR
        fq1_str=lambda w, input: ",".join(input.fq1),
        fq2_str=lambda w, input: ",".join(input.fq2) if input.fq2 and None not in input.fq2 else ""
    threads: 8
    conda:
        "../env/star.yaml"
    resources:
        mem_mb=50000,
        runtime=400
    shell:
        """
        # 1. Create the temp directory for STAR prefixes
        mkdir -p {params.prefix}

        # 2. Run STAR with ENCODE parameters
        # --outReadsUnmapped Fastx: Required to generate the unmapped fastq files
        # --readFilesCommand zcat: Required for .gz input
        # --quantMode GeneCounts TranscriptomeSAM: Required for both count tables and RSEM BAM
        # Multiple input files are separated by commas
        
        STAR --runThreadN {threads} \
            --genomeDir {input.idx} \
            --readFilesIn {params.fq1_str} {params.fq2_str} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.prefix} \
            --sjdbGTFfile {input.gtf} \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts TranscriptomeSAM \
            --outReadsUnmapped Fastx \
            --outFilterMultimapNmax 20 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --outFilterType BySJout \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --alignIntronMax 1000000 \
            --outSAMunmapped Within \
            --outSAMattributes NH HI AS NM MD \
            2> {log}

        # 3. Move and Rename Outputs to match Snakemake expectations
        
        # Genomic BAM
        mv {params.prefix}Aligned.sortedByCoord.out.bam {output.aln}
        
        # Transcriptome BAM
        mv {params.prefix}Aligned.toTranscriptome.out.bam {output.transcriptome_bam}
        
        # Gene Counts (ReadsPerGene.out.tab)
        mv {params.prefix}ReadsPerGene.out.tab {output.reads_per_gene}
        
        # Splice Junctions
        mv {params.prefix}SJ.out.tab {output.sj}
        
        # Log file
        mv {params.prefix}Log.final.out {output.log_final}

        # 4. Handle Unmapped files (Compress and Rename)
        # STAR outputs them as Unmapped.out.mate1 and Unmapped.out.mate2
        gzip -c {params.prefix}Unmapped.out.mate1 > {output.unmapped[0]}
        gzip -c {params.prefix}Unmapped.out.mate2 > {output.unmapped[1]}

        # 5. Cleanup temp directory
        rm -rf {params.prefix}
        """

# Note: merge_sample_runs rule is no longer needed since star_align
# now processes all runs together and outputs per-sample BAM files directly

rule samtools_index_bam:
    input:
        os.path.join(result_path, "Important_processed", "Bam", "{sample}_sortedByCoord.out.bam")
    output:
        os.path.join(result_path, "Important_processed", "Bam", "{sample}_sortedByCoord.out.bam.bai")
    log:
        os.path.join(result_path, "logs", "samtools", "index", "{sample}.log")
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
