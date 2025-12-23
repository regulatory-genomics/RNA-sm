import os

# Note: merge_transcriptome_bam rule is no longer needed since star_align
# now processes all runs together and outputs per-sample transcriptome BAM files directly

rule rsem_quant:
    """
    Quantify gene expression using RSEM.
    Matches ENCODE settings: Fixed seed, CI calculation, BAM input.
    """
    input:
        # CRITICAL: Must use the 'toTranscriptome' BAM from STAR (merged per sample)
        bam=os.path.join(result_path, "Important_processed", "Bam", "{sample}_Aligned.toTranscriptome.out.bam"),
        # RSEM Index marker file
        index_marker=os.path.join(config["resources"]["rsem_index_dir"], "rsem.seq")
    output:
        genes=os.path.join(result_path, "downstream_res", "rsem", "{sample}.genes.results"),
        isoforms=os.path.join(result_path, "downstream_res", "rsem", "{sample}.isoforms.results")
    params:
        seed=12345,
        fwd_prob=lambda w: get_rsem_forward_prob(w),
        # Detect paired-end from the sample info or assume paired
        paired = lambda w: get_rsem_paired_flag(w),
        index_prefix=os.path.join(config["resources"]["rsem_index_dir"], "rsem"),
        output_prefix=lambda w: os.path.join(result_path, "downstream_res", "rsem", w.sample)
    threads: 8
    resources:
        # RSEM needs significant RAM for CI calculation (approx 30GB+ for human)
        mem_mb=32000,
        runtime=120
    conda:
        "../env/rsem.yaml"
    shell:
        """
        rsem-calculate-expression --bam \
            --estimate-rspd \
            --calc-ci \
            --seed {params.seed} \
            -p {threads} \
            --no-bam-output \
            --forward-prob {params.fwd_prob} \
            {params.paired} \
            {input.bam} \
            {params.index_prefix} \
            {params.output_prefix}
        """

rule rsem_index:
    """
    Builds the RSEM index required for quantification.
    """
    input:
        fasta=config["resources"]["fasta"],  # e.g. genome.fa
        gtf=config["resources"]["gtf"]       # e.g. annotation.gtf
    output:
        # RSEM creates multiple files (rsem.grp, rsem.ti, rsem.seq, etc.)
        # We track the main configuration file to mark completion.
        seq=os.path.join(config["resources"]["rsem_index_dir"], "rsem.seq"),
        grp=os.path.join(config["resources"]["rsem_index_dir"], "rsem.grp")
    params:
        # The prefix includes the directory path + "rsem"
        prefix=os.path.join(config["resources"]["rsem_index_dir"], "rsem")
    threads: 8
    resources:
        mem_mb=16000
    conda:
        "../env/rsem.yaml"
    shell:
        """
        # Create the directory first
        mkdir -p $(dirname {params.prefix})
        
        # Build the index
        rsem-prepare-reference \
            --gtf {input.gtf} \
            --num-threads {threads} \
            {input.fasta} \
            {params.prefix}
        """