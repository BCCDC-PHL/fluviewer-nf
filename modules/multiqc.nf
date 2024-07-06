process multiqc {

    publishDir "${params.outdir}", pattern: "*_multiqc_report.html", mode:'copy'

    input:
    path '*'

    output:
    path '*_multiqc_*'

    script:
    """
    multiqc . -n ${params.run_name}_fluviewer-nf_multiqc_report.html
    """
}
