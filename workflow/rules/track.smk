import os


rule star_bedgraph:
    """
    Step 1: Generate bedGraph files from BAM using STAR in bedGraph mode.
    
    STAR command:
    --runMode inputAlignmentsFromBAM: Reads an existing BAM file (doesn't align)
    --outWigType bedGraph: Outputs bedGraph format
    --outWigStrand: Handles strandedness (Stranded or Unstranded)
    --outWigReferencesPrefix chr: Adds "chr" prefix to chromosome names
    
    Outputs (in directory):
    - Signal.Unique.str1.out.bg - Unique reads, strand 1 (always present)
    - Signal.UniqueMultiple.str1.out.bg - All reads, strand 1 (always present)
    - Signal.Unique.str2.out.bg - Unique reads, strand 2 (only for Stranded)
    - Signal.UniqueMultiple.str2.out.bg - All reads, strand 2 (only for Stranded)
    
    Note: Using directory output because unstranded libraries don't produce str2 files.
    """
    input:
        bam=lambda w: get_bam_for_tracks(w.sample),
        idx=config["resources"]["star_index"]
    output:
        # Use directory output instead of specific files
        # This handles the case where unstranded libraries don't produce str2 files
        outdir=directory(os.path.join(result_path, "middle_file", "bedgraph", "{sample}")),
        done=touch(os.path.join(result_path, "middle_file", "bedgraph", "{sample}", "star.done"))
    params:
        strandedness=lambda w: get_track_strandedness(w.sample),
    log:
        os.path.join(result_path, "logs", "star", "bedgraph", "{sample}.log")
    threads: 4
    resources:
        mem_mb=32000,
        runtime=120
    conda:
        "../env/star.yaml"
    shell:
        """
        mkdir -p {output.outdir}

        STAR --runMode inputAlignmentsFromBAM \
            --inputBAMfile {input.bam} \
            --genomeDir {input.idx} \
            --outWigType bedGraph \
            --outWigStrand {params.strandedness} \
            --outWigReferencesPrefix chr \
            --outFileNamePrefix {output.outdir}/ \
            2> {log}
        """


rule bedgraph_to_bigwig:
    """
    Step 2: Convert bedGraph to BigWig.
    """
    input:
        # 1. Connect strictly to the DIRECTORY output of star_bedgraph
        done_flag = os.path.join(result_path, "middle_file", "bedgraph", "{sample}", "star.done"),
        chrom_sizes = config["resources"]["chrom_sizes"]
    output:
        bw = os.path.join(result_path, "Important_processed", "Track", "{sample}_{suffix}.bw")
    params:
        # 2. Determine which file inside that directory we want (e.g. Signal.Unique.str1.out.bg)
        bg_dir = lambda w: os.path.join(result_path, "middle_file", "bedgraph", w.sample),
        filename = lambda w: get_star_bg_filename(w)
    log:
        os.path.join(result_path, "logs", "tracks", "{sample}_{suffix}.log")
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 60
    conda:
        "../env/kentutils.yaml"
    shell:
        """
        # Construct the full path here: {input.bg_dir}/{params.filename}
        
        bedSort {params.bg_dir}/{params.filename} {params.bg_dir}/{params.filename}.sorted 2> {log} && \
        
        bedGraphToBigWig {params.bg_dir}/{params.filename}.sorted {input.chrom_sizes} {output.bw} 2>> {log} && \
        
        rm {params.bg_dir}/{params.filename}.sorted
        """