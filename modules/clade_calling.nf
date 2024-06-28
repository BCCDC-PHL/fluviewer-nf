process clade_calling {

    conda "${projectDir}/environments/nextclade.yml"

    errorStrategy 'ignore'

    tag { sample_id }

    publishDir "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}/clade-calls", pattern: "${sample_id}_nextclade.*", mode:'copy'
    publishDir "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}/clade-calls", pattern: "${sample_id}*translation.fasta.gz", mode:'copy'

    input:
    tuple val(sample_id), path(consensus_seqs)

    output:
    tuple val(sample_id), path("*nextclade*"), emit: nextclade, optional: true
    tuple val(sample_id), path("${sample_id}_clade_calling_provenance.yml"), emit: provenance, optional: true

    script:
    """
    printf -- "process_name: nextclade\\n"  >> ${sample_id}_clade_calling_provenance.yml
    printf -- "tools:\\n"                   >> ${sample_id}_clade_calling_provenance.yml
    printf -- "  - tool_name: nextclade\\n" >> ${sample_id}_clade_calling_provenance.yml
    printf -- "    tool_version: \$(nextclade --version 2>&1  | cut -d ' ' -f 2)\\n" >> ${sample_id}_clade_calling_provenance.yml
    printf -- "    subcommand: run\\n"       >> ${sample_id}_clade_calling_provenance.yml

    [ ! -f  ${sample_id}_HA_consensus.fa ] && ln -sf *HA_consensus.fa ${sample_id}_HA_consensus.fa

    FOUND=true

    if [ `grep "H1" ${sample_id}_HA_consensus.fa` ]; then
        dataset=${params.h1_dataset}
    elif [ `grep "H3" ${sample_id}_HA_consensus.fa` ]; then 
        dataset=${params.h3_dataset}
    elif [ `grep "H5" ${sample_id}_HA_consensus.fa` ]; then 
        dataset=${params.h5_dataset}
    else 
        echo "WARNING: None of H1, H3, or H5 were detected in the HA consensus file. No dataset available. Exiting."
        FOUND=false
        dataset="NONE"
    fi 

    if [ \$FOUND == true ]; then 

        nextclade run --input-dataset \$dataset \
        --output-fasta=${sample_id}_nextclade.aligned.fasta.gz \
        --output-json=${sample_id}_nextclade.json \
        --output-ndjson=${sample_id}_nextclade.ndjson \
        --output-csv=${sample_id}_nextclade.csv \
        --output-tsv=${sample_id}_nextclade.tsv \
        --output-tree=${sample_id}_nextclade_auspice.json \
        --output-translations=${sample_id}_nextclade_{gene}.translation.fasta.gz \
        ${sample_id}_HA_consensus.fa

        LOCATION="\${dataset}\\n" 
        VERSION="\$(grep "tag" \${dataset}/tag.json)\\n" 
    else
        LOCATION="NONE_INVALID_HA_TYPE\\n"
        VERSION="NONE_INVALID_HA_TYPE\\n"
    fi 

    """
}
