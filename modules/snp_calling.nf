process snp_calling {

    errorStrategy 'ignore'

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/snp-calls", pattern: "*mutations.tsv", mode:'copy'
    publishDir "${params.outdir}/${sample_id}/snp-calls/pairwise", pattern: "*fa", mode:'copy'

    input:
    tuple val(sample_id), path(consensus_seqs)
    tuple path(blastx_db_path), val(blastx_db_name)

    output:
    tuple val(sample_id), path("${sample_id}_*mutations.tsv"), emit: mutations
    tuple val(sample_id), path("${sample_id}_snp_calling_provenance.yml"), emit: provenance, optional: true

    script:
    """
    printf -- "- process_name: snp_calling\\n" >> ${sample_id}_snp_calling_provenance.yml
    printf -- "  tools:\\n"                   >> ${sample_id}_snp_calling_provenance.yml
    printf -- "    - tool_name: blastx\\n"    >> ${sample_id}_snp_calling_provenance.yml
    printf -- "      tool_version: \$(blastx -version | head -n 1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_snp_calling_provenance.yml
    printf -- "  databases:\\n"               >> ${sample_id}_snp_calling_provenance.yml
    printf -- "    - database_name: ${blastx_db_name}\\n" >> ${sample_id}_snp_calling_provenance.yml

    aa_mutation_calling.py \
    --database ${blastx_db_path}/${blastx_db_name} \
    --query ${consensus_seqs} \
    --output ${sample_id}_{segment}_mutations.tsv \
    --outaln ${sample_id}_{segment}_aln.fa
    """
}
