process genoflu {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/", pattern: "*tsv", mode: 'copy'


    input:
    tuple val(sample_id), path(consensus_seqs)

    output:
    tuple val(sample_id), path("${sample_id}_genoflu.tsv"), emit: tsv
    tuple val(sample_id), path("${sample_id}*_provenance.yml"), emit: provenance

    script:
    """
    GENOFLU_VERSION=\$(genoflu.py --version | cut -d' ' -f3)
    GENOFLU_DB_PATH="\$(which genoflu.py | xargs dirname | xargs dirname)/dependencies"

    printf -- "- process_name: genoflu\\n" >> ${sample_id}_genoflu_provenance.yml
    printf -- "  tools:\\n"                >> ${sample_id}_genoflu_provenance.yml
    printf -- "    - tool_name: genoflu\\n"    >> ${sample_id}_genoflu_provenance.yml
    printf -- "      tool_version: \$GENOFLU_VERSION\\n" >> ${sample_id}_genoflu_provenance.yml

    genoflu.py \
	-f ${consensus_seqs} \
	-i \${GENOFLU_DB_PATH}/fastas/ \
	-c \${GENOFLU_DB_PATH}/genotype_key.xlsx \
	-n ${sample_id}

    mv ${sample_id}*tsv ${sample_id}_genoflu_tmp.txt

    echo -e "Genoflu Version\n\${GENOFLU_VERSION}" > gf_version_column.txt

    paste ${sample_id}_genoflu_tmp.txt gf_version_column.txt > ${sample_id}_genoflu.tsv
    
    """
}
