process snp_calling {

    errorStrategy 'ignore'

    tag { sample_id }

    publishDir "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}/snp-calls", pattern: "*mutations.tsv", mode:'copy'
    publishDir "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}/snp-calls/pairwise", pattern: "*fa", mode:'copy'

    input:
    tuple val(sample_id), path(consensus_seqs)
    tuple path(blastx_db_path), val(blastx_db_name)

    output:
    tuple val(sample_id), path("${sample_id}_*mutations.tsv"), emit: mutations
    tuple val(sample_id), path("${sample_id}/*blastx.tsv"), emit: blast, optional: true
    // tuple val(sample_id), path("${sample_id}_clade_provenance.yml"), emit: provenance, optional: true

    script:
    """ 
    export BLASTDB="${blastx_db_path}"
    blastx -query ${consensus_seqs} -db ${blastx_db_name} -outfmt 6 > ${sample_id}_blastx.tsv && 

    snp-calling.py \
    --reference_db ${blastx_db_path}/${blastx_db_name} \
    --blast ${sample_id}_blastx.tsv \
    --query ${consensus_seqs} \
    --output ${sample_id}_{segment}_mutations.tsv \
    --outaln ${sample_id}_{segment}_aln.fa
    """
}
