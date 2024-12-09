process make_coverage_plot {

    errorStrategy 'ignore'

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/", pattern: "*png", mode:'copy'

    input:
    tuple val(sample_id), path(samtools_depth)

    output:
    tuple val(sample_id), path("${sample_id}_*png"), emit: plot

    script:
    """
    coverage_plot.py \
    --input ${samtools_depth} \
	--output ${sample_id}_coverage.png
    """
}