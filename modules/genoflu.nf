process genoflu {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/", pattern: "*tsv", mode: 'copy'


    input:
    tuple val(sample_id), path(consensus_seqs), path(genoflu_db)

    output:
    tuple val(sample_id), path("${sample_id}_genoflu.tsv"), emit: tsv
    tuple val(sample_id), path("${sample_id}*_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: genoflu\\n" >> ${sample_id}_genoflu_provenance.yml
    printf -- "  tools:\\n"                >> ${sample_id}_genoflu_provenance.yml
    printf -- "    - tool_name: genoflu\\n"    >> ${sample_id}_genoflu_provenance.yml
    printf -- "      tool_version: \$(genoflu.py --version | cut -d' ' -f3)\\n" >> ${sample_id}_genoflu_provenance.yml

    ${genoflu_db}/bin/genoflu.py \
	-f ${consensus_seqs} \
	-i ${genoflu_db}/dependencies/fastas/ \
	-c ${genoflu_db}/dependencies/genotype_key.xlsx \
	-n ${sample_id}

    mv ${sample_id}*tsv ${sample_id}_genoflu.tsv
    """
}
