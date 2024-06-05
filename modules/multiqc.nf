process multiqc {

    publishDir "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}", pattern: "*_multiqc_report.html", mode:'copy'

    input:
    path '*'

    output:
    path '*_multiqc_*'

    script:
    """
    multiqc . -n ${params.run_name}_multiqc_report.html
    """
}
