import os

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
        json = os.path.join(result_path, "Report", "gene_type_counts", "{sample}.json")
    params:
        # Path to the python script you saved
        script = "workflow/scripts/rna_qc.py"
    resources:
        mem_mb = 4000,   # Parsing BAMs requires decent memory, but less than RSEM
        runtime = 60     # Should be fast (< 30 mins usually)
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
            > {log} 2>&1
        """

rule multiqc:
    input:
        expand(
            os.path.join(result_path, "Report", "fastp", "{sample_run}_fastp.json"),
            sample_run=annot.index
        ),
        expand(
            os.path.join(result_path, "Report", "star", "{sample_run}_Log.final.out"),
            sample_run=annot.index
        ),
    output:
        html=os.path.join(result_path, "Report", "multiqc_report.html"),
        data=directory(os.path.join(result_path, "Report", "multiqc_report_data"))
    log:
        os.path.join(result_path, "logs", "multiqc.log")
    conda:
        "../env/multiqc.yaml"
    params:
        outdir=os.path.join(result_path, "Report")
    shell:
        "multiqc {input} -o {params.outdir} --filename multiqc_report.html --force 2> {log}"
