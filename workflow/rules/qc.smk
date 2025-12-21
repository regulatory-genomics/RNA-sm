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
