import os

rule fastp:
    input:
        r1 = lambda w: get_units_fastqs(w)[0],
        r2 = lambda w: get_units_fastqs(w)[1],
    output:
        r1 = temp(os.path.join(result_path, "trimmed", "{sample_run}_1.fq.gz")),
        r2 = temp(os.path.join(result_path, "trimmed", "{sample_run}_2.fq.gz")),
        report_html = os.path.join(result_path, "qc", "fastp", "{sample_run}_fastp.html"),
        report_json = os.path.join(result_path, "qc", "fastp", "{sample_run}_fastp.json"),
    conda: 
        "../env/fastp.yaml"
    resources:
        mem_mb=16000,
        runtime = 60,
    log:
        "logs/fastp/{sample_run}.log"
    threads: 4
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
        "--detect_adapter_for_pe --trim_poly_g --thread {threads} "
        "-j {output.report_json} -h {output.report_html} 2> {log}"

rule multiqc:
    input:
        expand(
            os.path.join(result_path, "qc", "fastp", "{sample_run}_fastp.json"),
            sample_run=annot.index
        ),
        expand(
            os.path.join(result_path, "align", "{sample_run}_Log.final.out"),
            sample_run=annot.index
        ),
    output:
        html=os.path.join(result_path, "qc", "multiqc_report.html"),
        data=directory(os.path.join(result_path, "qc", "multiqc_report_data"))
    log:
        "logs/multiqc.log"
    conda:
        "../env/multiqc.yaml"
    params:
        outdir=os.path.join(result_path, "qc")
    shell:
        "multiqc {input} -o {params.outdir} --filename multiqc_report.html --force 2> {log}"
