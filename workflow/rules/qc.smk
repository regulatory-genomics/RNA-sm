rule fastqc:
    input:
        fq1=lambda w: get_units_fastqs(w)[0],
        fq2=lambda w: get_units_fastqs(w)[1]
    output:
        # Name the outputs to match the inputs
        html1="results/fastqc/{sample}_{unit}_1_fastqc.html",
        html2="results/fastqc/{sample}_{unit}_2_fastqc.html",
        zip1="results/fastqc/{sample}_{unit}_1_fastqc.zip",
        zip2="results/fastqc/{sample}_{unit}_2_fastqc.zip"
    params:
        outdir="results/fastqc"
    log:
        "logs/fastqc/{sample}-{unit}.log"
    shell:
                """
        # 1. Run fastqc on both input files
        fastqc {input.fq1} {input.fq2} --outdir {params.outdir} &> {log}

        # 2. Figure out the default output names fastqc created
        # The 'basename' command strips directory paths and optionally a suffix
        fq1_base=$(basename {input.fq1} .fq.gz)
        fq2_base=$(basename {input.fq2} .fq.gz)

        # 3. Rename the actual files to the names Snakemake expects
        mv {params.outdir}/${{fq1_base}}_fastqc.html {output.html1}
        mv {params.outdir}/${{fq2_base}}_fastqc.html {output.html2}
        mv {params.outdir}/${{fq1_base}}_fastqc.zip {output.zip1}
        mv {params.outdir}/${{fq2_base}}_fastqc.zip {output.zip2}
        """


rule trim_galore_pe:
    input:
        get_units_fastqs,
    output:
        fasta_fwd="results/trimmed/{sample}-{unit}_1.fq.gz",
        report_fwd="results/trimmed/reports/{sample}-{unit}_1_trimming_report.txt",
        fasta_rev="results/trimmed/{sample}-{unit}_2.fq.gz",
        report_rev="results/trimmed/reports/{sample}-{unit}_2_trimming_report.txt",
    threads: 1
    log:
        "logs/trim_galore/{sample}-{unit}.log"
    wrapper:
        "v7.6.0/bio/trim_galore/pe"
