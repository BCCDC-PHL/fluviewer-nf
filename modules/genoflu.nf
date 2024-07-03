process genoflu {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/", pattern: "*tsv", mode: 'copy'


    input:
    tuple val(sample_id), path(consensus_seqs), path(genoflu_path)

    output:
    tuple val(sample_id), path("${sample_id}_genoflu.tsv"), emit: tsv
    tuple val(sample_id), path("${sample_id}*_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: genoflu\\n" >> ${sample_id}_genoflu_provenance.yml
    printf -- "  tools:\\n"                >> ${sample_id}_genoflu_provenance.yml
    printf -- "    - tool_name: genoflu\\n"    >> ${sample_id}_genoflu_provenance.yml
    printf -- "      tool_version: \$(genoflu.py --version | cut -d' ' -f3)\\n" >> ${sample_id}_genoflu_provenance.yml

    genoflu.py \
	-f ${consensus_seqs} \
	-i ${genoflu_path}/dependencies/fastas/ \
	-c ${genoflu_path}/dependencies/genotype_key.xlsx \
	-n ${sample_id}

    mv ${sample_id}*tsv ${sample_id}_genoflu.tsv
    """
}

process pull_genoflu {
    
    executor 'local'
    storeDir "${params.genoflu_cache}"

    input:
    val(genoflu_github_url)

    output:
    path("GenoFLU"), emit: repo

    script:
    """
    git clone ${genoflu_github_url} GenoFLU
    """
}

process checkout_genoflu {

    input:
    path(genoflu_path)
    val(genoflu_version)

    output:
    path("checkout_log.txt"), emit: log
    val("\${VERSION}"), emit: version

    script:
    """
    CWD=\$PWD
    cd ${genoflu_path}

    echo "Requested GenoFLU Version: ${genoflu_version}" > \${CWD}/checkout_log.txt

    git pull origin main

    LATEST_VERSION=`git describe --tags \$(git rev-list --tags --max-count=1)`

    if [ ${genoflu_version} == "LATEST" ]; then
        VERSION=\${LATEST_VERSION}
    else
        VERSION=${genoflu_version}
    fi 

    echo "Translated GenoFLU Version: \${VERSION}" >> \${CWD}/checkout_log.txt
    
    TAG_LIST=(\$(git tag))

    echo "Available Tags: \${TAG_LIST[@]}" >> \${CWD}/checkout_log.txt

    if ! [[ \${TAG_LIST[@]} =~ (^|[[:space:]])\$VERSION(\$|[[:space:]]) ]]; then 
        echo "ERROR: The provided version does not exist on the GitHub repo after pulling. Exiting."
        exit 1
    fi 

    git checkout tags/\${VERSION} 2>> \${CWD}/checkout_log.txt

    """
}
