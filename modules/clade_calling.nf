process clade_calling {

    conda "${projectDir}/environments/nextclade.yml"

    errorStrategy 'ignore'

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/clade-calls", pattern: "${sample_id}_nextclade.*", mode:'copy'
    publishDir "${params.outdir}/${sample_id}/clade-calls", pattern: "${sample_id}*translation.fasta.gz", mode:'copy'

    input:
    tuple val(sample_id), path(ha_consensus_seq)

    output:
    tuple val(sample_id), path("*nextclade*"), emit: nextclade
    tuple val(sample_id), path("${sample_id}_clade_calling_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: nextclade\\n"  >> ${sample_id}_clade_calling_provenance.yml
    printf -- "  tools:\\n"                   >> ${sample_id}_clade_calling_provenance.yml
    printf -- "    - tool_name: nextclade\\n" >> ${sample_id}_clade_calling_provenance.yml
    printf -- "      tool_version: \$(nextclade --version 2>&1  | cut -d ' ' -f 2)\\n" >> ${sample_id}_clade_calling_provenance.yml
    printf -- "      subcommand: run\\n"       >> ${sample_id}_clade_calling_provenance.yml


    if [ `grep "H1" ${ha_consensus_seq}` ]; then
        if [ ${params.h1_dataset} == "NO_FILE" ]; then
            echo "WARNING: H1 subtype detected in the HA consensus file, but no H1 dataset provided. Please provide an H1 nextclade dataset with the --h1_dataset flag."
            exit 10
        else
            dataset=${params.h1_dataset}
        fi
    elif [ `grep "H3" ${ha_consensus_seq}` ]; then 
        if [ ${params.h3_dataset} == "NO_FILE" ]; then
            echo "WARNING: H3 subtype detected in the HA consensus file, but no H3 dataset provided. Please provide an H3 nextclade dataset with the --h3_dataset flag."
            exit 10
        else
            dataset=${params.h3_dataset}
        fi
    elif [ `grep "H5" ${ha_consensus_seq}` ]; then
        if [ ${params.h5_dataset} == "NO_FILE" ]; then
            echo "WARNING: H5 subtype detected in the HA consensus file, but no H5 dataset provided. Please provide an H5 nextclade dataset with the --h5_dataset flag."
            exit 10
        else
            dataset=${params.h5_dataset}
        fi
    else 
        echo "WARNING: None of H1, H3, or H5 were detected in the HA consensus file. Exiting."
        exit 10
    fi

    printf -- "  databases:\\n"                >> ${sample_id}_clade_calling_provenance.yml
    printf -- "    - database_name: \$(grep 'dbname' \${dataset}/metadata.json | cut -d'\"' -f4) \\n"                  >> ${sample_id}_clade_calling_provenance.yml
    printf -- "      database_path: \${dataset} \\n"                                                      >> ${sample_id}_clade_calling_provenance.yml
    printf -- "      database_sha256: \$(sha256sum \${dataset}/tree.json | cut -d' ' -f1) \\n"                >> ${sample_id}_clade_calling_provenance.yml


    nextclade run \
        --input-dataset \$dataset \
        --output-fasta=${sample_id}_nextclade.aligned.fasta.gz \
        --output-json=${sample_id}_nextclade.json \
        --output-ndjson=${sample_id}_nextclade.ndjson \
        --output-csv=${sample_id}_nextclade.csv \
        --output-tsv=${sample_id}_nextclade.tsv \
        --output-tree=${sample_id}_nextclade_auspice.json \
        --output-translations=${sample_id}_nextclade_{cds}.translation.fasta.gz \
        ${ha_consensus_seq}
    """
}
