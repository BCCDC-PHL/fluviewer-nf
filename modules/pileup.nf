process make_pileup {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "*pileup.tsv",  mode:'copy'

    input:
    tuple val(sample_id), path(bam), path(bam_index)

    output:
    tuple val(sample_id), path("${sample_id}_pileup.tsv"), emit: pileup

    script:
    """
    make_pileup.py --input ${bam} --output ${sample_id}_pileup.tsv
    """
}

process plot_pileup {

    tag { sample_id }

    conda "${projectDir}/environments/fluviewer.yml"

    publishDir "${params.outdir}/${sample_id}/plots", pattern: "*png",  mode:'copy'

    input:
    tuple val(sample_id), path(pileup)

    output:
    tuple val(sample_id), path("${sample_id}_*combined.png"), emit: combined
    tuple val(sample_id), path("${sample_id}_*{NA,HA,PB2,PB1,PA,NS,M,NP}.png"), emit: separate

    script:
    """
    plot_pileup.py --input ${pileup} --output ${sample_id}_pileup_combined.png -m facet

    plot_pileup.py --input ${pileup} --output ${sample_id}_pileup.png -m separate
    """
}