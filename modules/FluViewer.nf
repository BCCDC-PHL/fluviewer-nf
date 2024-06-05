process fluviewer {

    tag { sample_id }

    memory  { 50.GB * task.attempt }
    errorStrategy { (task.exitStatus == 2 && task.attempt <= maxRetries) ? 'retry' : 'ignore' } 
    maxRetries 5

    publishDir "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}", pattern: "${sample_id}_fluviewer/${sample_id}*", mode:'copy', saveAs: { filename -> filename.split("/").last() }
    publishDir "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}", pattern: "${sample_id}_fluviewer/*tsv", mode:'copy', saveAs: { filename -> filename.split("/").last() }
    publishDir "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}", pattern: "${sample_id}_fluviewer/spades_output", mode:'copy', saveAs: { filename -> "spades_output" }
    publishDir "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}", pattern: ".*", mode:'copy'
    publishDir "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}", pattern: "${sample_id}_fluviewer/logs", mode:'copy', saveAs: { filename -> "fluviewer_logs" }
    publishDir "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}", pattern: ".exitcode", mode:'copy'
    publishDir "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}", pattern: ".command.*", mode:'copy'
  
    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(db)

    output:
    tuple val(sample_id), path("${sample_id}_fluviewer/${sample_id}*.bam"), emit: alignment
    tuple val(sample_id), path("${sample_id}_fluviewer/${sample_id}*.bam.bai"), emit: alignmentindex, optional: true
    tuple val(sample_id), path("${sample_id}_fluviewer/${sample_id}*report.tsv"), emit: reports, optional: true
    tuple val(sample_id), path("${sample_id}_fluviewer/${sample_id}*_consensus.fa"), emit: consensus_seqs, optional: true
    tuple val(sample_id), path("${sample_id}_fluviewer/${sample_id}*consensus_seqs.fa"), emit: consensus_main
    tuple val(sample_id), path("${sample_id}_fluviewer/${sample_id}*_HPAI.tsv"), emit: HPAI, optional: true
    tuple val(sample_id), path("${sample_id}_fluviewer/${sample_id}*_cov.png"), emit: coverage_plot, optional: true
    tuple val(sample_id), path("${sample_id}_fluviewer/${sample_id}*_variants.vcf"), emit: vcf, optional: true
    tuple val(sample_id), path("${sample_id}_fluviewer/logs"), emit: fluviewer_logs
    tuple val(sample_id), path("${sample_id}_FluViewer_provenance.yml"), emit: provenance
    tuple val(sample_id), path("${sample_id}_fluviewer/${sample_id}*_mapping_refs.fa"), emit: ref_seqs_for_mapping, optional: true
    tuple val(sample_id), path("${sample_id}_fluviewer/contigs_blast.tsv"), emit: contig_blast_results, optional: true
    tuple val(sample_id), path("${sample_id}_fluviewer/spades_output"), emit: spades_results, optional: true
    tuple val(sample_id), path("${sample_id}_fluviewer/${sample_id}*.png"), emit: depth_cov_plot, optional: true

    script:
    garbage_collection = params.keep_interfiles ? '-g' : ''
    """
    printf -- "- process_name: fluviewer\\n"                  >> ${sample_id}_FluViewer_provenance.yml
    printf -- "  tools:\\n"                                   >> ${sample_id}_FluViewer_provenance.yml
    printf -- "    - tool_name: FluViewer\\n"                 >> ${sample_id}_FluViewer_provenance.yml
    printf -- "      tool_version: \$(FluViewer | sed -n '4 p')\\n" >> ${sample_id}_FluViewer_provenance.yml
    printf -- "  databases:\\n"                               >> ${sample_id}_FluViewer_provenance.yml
    printf -- "    - database_name: ${db}\\n"                 >> ${sample_id}_FluViewer_provenance.yml
    printf -- "      database_path: \$(readlink -f ${db})\\n" >> ${sample_id}_FluViewer_provenance.yml
    printf -- "      database_sha256: \$(shasum -a 256 ${db}|awk '{print \$1}')\\n" >> ${sample_id}_FluViewer_provenance.yml
  
    EXITCODE=0
    (FluViewer \
	${garbage_collection} \
	-T ${task.cpus} \
	-f ${reads_1} -r ${reads_2} \
	-n ${sample_id}_fluviewer \
	-d ${db} \
	-D ${params.min_depth}  \
	-q ${params.min_q} \
	-i ${params.min_ident} \
	-M 40 && EXITCODE=\$?) \
	|| EXITCODE=\$?

    echo "Extracting NA and HA consensus sequences..."

    if [ `grep "|HA|" ${sample_id}_fluviewer/${sample_id}*consensus_seqs.fa` ]; then 
        grep -A1 "|HA|" ${sample_id}_fluviewer/${sample_id}*consensus_seqs.fa > ${sample_id}_fluviewer/${sample_id}_HA_consensus.fa
    else
        echo "No HA consensus sequence generated."
    fi

    if [ `grep "|NA|" ${sample_id}_fluviewer/${sample_id}*consensus_seqs.fa` ]; then 
        grep -A1 "|NA|" ${sample_id}_fluviewer/${sample_id}*consensus_seqs.fa > ${sample_id}_fluviewer/${sample_id}_NA_consensus.fa
    else
        echo "No NA consensus sequence generated."
    fi

    if [[ ! -f ${sample_id}_fluviewer/${sample_id}_HA_consensus.fa ]]; then
        echo "HA segment consensus not generated. Skipping FindCleave.py..."
    else
	python ${projectDir}/bin/FindCleave.py -i ${sample_id}_fluviewer/${sample_id}_HA_consensus.fa -o ${sample_id}_fluviewer/${sample_id}_HPAI.tsv
        echo "Finished running FindCleave.py."
    fi

  echo \$EXITCODE > .exitcode

  OUTPATH=${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}
  if [[ \$OUTPATH != /* ]]; then 
    OUTPATH=${workflow.launchDir}/\$OUTPATH
  fi

  cp .command.* \$OUTPATH
  cp .exitcode \$OUTPATH
  exit \$EXITCODE
  """
}
