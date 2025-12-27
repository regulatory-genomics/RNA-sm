import os
import itertools
import json
import pandas as pd

rule fastp:
    input:
        r1 = lambda w: get_units_fastqs(w)[0],
        r2 = lambda w: get_units_fastqs(w)[1],
    output:
        r1 = temp(os.path.join(result_path, "middle_file", "trimmed", "{sample_run}_1.fq.gz")),
        r2 = temp(os.path.join(result_path, "middle_file", "trimmed", "{sample_run}_2.fq.gz")),
        report_html = os.path.join(result_path, "Report", "fastp", "{sample_run}_fastp.html"),
        report_json = os.path.join(result_path, "Report", "fastp", "{sample_run}_fastp.json"),
    conda: 
        "../env/fastp.yaml"
    resources:
        mem_mb=16000,
        runtime = 60,
    log:
        os.path.join(result_path, "logs", "fastp", "{sample_run}.log")
    threads: 4
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
        "--detect_adapter_for_pe --trim_poly_g --thread {threads} "
        "-j {output.report_json} -h {output.report_html} 2> {log}"

rule quantify_gene_expression:
    input:
        rsem_tsv = os.path.join(result_path, "downstream_res", "rsem", "{sample}.genes.results"),
    output:
        gene_exp_json = os.path.join(result_path, "Report", "rsem", "{sample}.genes.json"),
    conda:
        "../env/pandas.yaml"
    run:
        import pandas as pd
        def calculate_number_of_genes_detected(quant_tsv, threshold_of_detection=1):
            """
            Calculates number of rows where value on the column TPM is greater than
            the threshold_of_detection.
            Args:
                quant_tsv: filename of a .tsv of RSEM quants
            Returns:
                int number_of_genes_detected: number of entries > threshold_of_detection
                in TPM column
            """
            quants = pd.read_csv(quant_tsv, sep="\t")
            number_of_genes_detected = sum(quants["TPM"] > threshold_of_detection)
            return number_of_genes_detected

        num_genes = calculate_number_of_genes_detected(input.rsem_tsv)
        result = {"num_genes_detected": int(num_genes)}
        import json
        with open(output.gene_exp_json, "w") as f:
            json.dump(result, f, indent=2)

rule qc_gene_type_count:
    """
    Calculate reads per gene type (protein_coding, rRNA, etc.)
    using the ENCODE rna_qc.py script.
    """
    input:
        # Uses the exact same transcriptome BAM as RSEM
        bam = os.path.join(result_path, "Important_processed", "Bam", "{sample}_Aligned.toTranscriptome.out.bam"),
        # The TSV mapping file defined in config.yaml
        tsv = config["resources"]["gene_type_tsv"]
    output:
        # Output JSON file with counts
        json = os.path.join(result_path, "Report", "gene_type_counts", "{sample}_gene_type_count.json")
    params:
        # Path to the python script you saved
        script = "workflow/scripts/rna_qc.py"
    resources:
        mem_mb = 4000,
        runtime = 60     
    conda:
        "../env/pysam.yaml" # Must contain pysam
    log:
        os.path.join(result_path, "logs", "qc_gene_type", "{sample}.log")
    shell:
        """
        python {params.script} \
            --input_bam {input.bam} \
            --tr_id_to_gene_type_tsv {input.tsv} \
            --output_filename {output.json} \
            --log_file {log} \
            > {log} 2>&1
        """

rule rnaseqqc:
    """
    Run rnaseqc quality control tool on a BAM file using a GTF annotation.
    Produces multiple output files with QC metrics, gene and exon counts, TPM, etc.
    Handles optional BED file only if provided.
    The BAM file is the sorted BAM file from STAR alignment.
    """
    input:
        bam=os.path.join(result_path, "Important_processed", "Bam", "{sample}_sortedByCoord.out.bam"),
        gtf=config["resources"]["collapsed_gtf"],
    output:
        metrics_tsv=os.path.join(result_path, "Report", "rnaseqqc", "{sample}.metrics.tsv"),
        exon_reads_gct=os.path.join(result_path, "Report", "rnaseqqc", "{sample}.exon_reads.gct"),
        gene_reads_gct=os.path.join(result_path, "Report", "rnaseqqc", "{sample}.gene_reads.gct"),
        gene_tpm_gct=os.path.join(result_path, "Report", "rnaseqqc", "{sample}.gene_tpm.gct"),
        coverage_tsv=os.path.join(result_path, "Report", "rnaseqqc", "{sample}.coverage.tsv")
    params:
        outdir=os.path.join(result_path, "Report", "rnaseqqc"),
        sample="{sample}",
        bed_flag=(f"--bed {config['resources'].get('rnaseqc_bed')}"
                if config.get("resources", {}).get("rnaseqc_bed") else "")
    log:
        os.path.join(result_path, "logs", "rnaseqqc", "{sample}.log")
    conda:
        "../env/rnaseqqc.yaml"
    shell:
        r"""
        mkdir -p {params.outdir}
        rnaseqc \
            {input.gtf} \
            {input.bam} \
            {params.outdir} \
            --sample {params.sample} \
            {params.bed_flag} \
            --coverage \
            > {log} 2>&1
        """

rule collapse_gtf:
    input:
        gtf=config["resources"]["gtf"],
    output:
        collapsed_gtf=config["resources"]["collapsed_gtf"],
    conda:
        "../env/collapse_gtf.yaml"
    log:
        os.path.join(result_path, "logs", "collapse_gtf.log")
    shell:
        "python workflow/scripts/collapse_gene.py {input.gtf} {output.collapsed_gtf} > {log} 2>&1"

rule mad_qc:
    input:
        # Dynamically fetch ALL samples for this replicate_name/group
        files = lambda w: [
            os.path.join(result_path, "downstream_res", "rsem", f"{s}.genes.results")
            for s in get_replicate_samples(w.replicate_name)
        ]
    output:
        # Directory for individual pairwise files
        out_dir = directory(os.path.join(result_path, "Report", "mad_qc", "{replicate_name}")),
        # New: The final summary file with averaged stats
        summary = os.path.join(result_path, "Report", "mad_qc", "{replicate_name}.summary_metrics.json")
    conda:
        "../env/pandas.yaml"
    log:
        os.path.join(result_path, "logs", "mad_qc_{replicate_name}.log")
    run:
        # 1. Setup output directory
        if not os.path.exists(output.out_dir):
            os.makedirs(output.out_dir)
            
        samples = get_replicate_samples(wildcards.replicate_name)
        files_map = dict(zip(samples, input.files))
        
        # List to collect data frames or dicts for aggregation
        all_metrics = []
        
        # 2. Iterate over all unique pairs
        # itertools.combinations('ABCD', 2) --> AB AC AD BC BD CD
        for s1, s2 in itertools.combinations(samples, 2):
            
            file1 = files_map[s1]
            file2 = files_map[s2]
            pair_id = f"{s1}_vs_{s2}"
            
            # Define filenames for this specific pair
            json_out = os.path.join(output.out_dir, f"{pair_id}.metrics.json")
            plot_out = os.path.join(output.out_dir, f"{pair_id}.plot.png")
            
            # 3. Run the QC script for this pair
            shell(
                f"""
                python workflow/scripts/mad_qc.py {file1} {file2} \
                    --json {json_out} \
                    --plot {plot_out} \
                    2>> {log}
                """
            )
            
            # 4. Read the results back immediately
            try:
                with open(json_out, 'r') as f:
                    data = json.load(f)
                    # Optional: Add metadata so you know which pair this was
                    data['pair'] = pair_id 
                    all_metrics.append(data)
            except Exception as e:
                # Log error if a specific pair fails, but don't crash the whole pipeline immediately
                log_path = log if isinstance(log, str) else log[0]
                with open(log_path, 'a') as l:
                    l.write(f"\n[Error] Could not read JSON for {pair_id}: {e}\n")
        
        # 5. Calculate Means and Save Summary
        if all_metrics:
            # Create a Pandas DataFrame
            df = pd.DataFrame(all_metrics)
            
            # Calculate mean of numeric columns (automatically ignores 'pair' string column)
            summary_stats = df.mean(numeric_only=True).to_dict()
            
            # Add a count of how many pairs were used
            summary_stats['num_pairs_evaluated'] = len(all_metrics)
            
            # Save the aggregated summary
            with open(output.summary, 'w') as f:
                json.dump(summary_stats, f, indent=4)
        else:
            # Handle case with single replicate (no pairs possible)
            with open(output.summary, 'w') as f:
                json.dump({"status": "No pairs available for aggregation"}, f, indent=4)

rule multiqc:
    input:
        # Keep specific files for dependency tracking
        expand(
            os.path.join(result_path, "Report", "fastp", "{sample_run}_fastp.json"),
            sample_run=annot.index
        ),
        expand(
            os.path.join(result_path, "Report", "star", "{sample}_Log.final.out"),
            sample=samples.keys()
        ),
        # Add plugin input files for dependency tracking
        expand(
            os.path.join(result_path, "Report", "rnaseqqc", "{sample}.metrics.tsv"),
            sample=samples.keys()
        ),
        expand(
            os.path.join(result_path, "Report", "gene_type_counts", "{sample}_gene_type_count.json"),
            sample=samples.keys()
        ),
        expand(
            os.path.join(result_path, "Report", "rsem", "{sample}.genes.json"),
            sample=samples.keys()
        ),
        expand(
            os.path.join(result_path, "Report", "mad_qc", "{replicate_name}"),
            replicate_name=replicate_samples
        ),
        # Include sample CSV as input so it's available for the plugin
        sample_csv=config["samples"]
    output:
        html=os.path.join(result_path, "Report", "multiqc_report.html"),
        data=directory(os.path.join(result_path, "Report", "multiqc_report_data"))
    log:
        os.path.join(result_path, "logs", "multiqc.log")
    conda:
        "../env/multiqc.yaml"
    params:
        outdir=os.path.join(result_path, "Report"),
        report_dir=os.path.join(result_path, "Report"),
        # Resolve absolute path to sample CSV for environment variable
        sample_csv_abs=os.path.abspath(config["samples"])
    shell:
        # Copy sample CSV to Report directory and set environment variable for plugin
        # Search the entire Report directory so MultiQC plugin can find all files
        r"""
        # Copy sample CSV to Report directory as backup location
        cp {input.sample_csv} {params.outdir}/sample_annotation.csv
        
        # Set environment variable for MultiQC plugin
        export MULTIQC_SAMPLE_SHEET={params.sample_csv_abs}
        
        # Run MultiQC
        multiqc {params.report_dir} -o {params.outdir} --filename multiqc_report.html --force 2> {log}
        """
