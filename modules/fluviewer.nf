process fluviewer {

    tag { sample_id }

    memory  { 50.GB * task.attempt }
    errorStrategy { (task.exitStatus == 2 && task.attempt <= maxRetries) ? 'retry' : 'ignore' } 
    maxRetries 5

    publishDir "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}", pattern: "${sample_id}*", mode:'copy', saveAs: { filename -> filename.split("/").last() }
    publishDir "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}", pattern: "*tsv", mode:'copy', saveAs: { filename -> filename.split("/").last() }
    //publishDir "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}", pattern: "${sample_id}_fluviewer/spades_output", mode:'copy', saveAs: { filename -> "spades_output" }
    publishDir "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}", pattern: ".*", mode:'copy'
    publishDir "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}", pattern: "logs", mode:'copy', saveAs: { filename -> "fluviewer_logs" }
    publishDir "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}", pattern: ".exitcode", mode:'copy'
    publishDir "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}", pattern: ".command.*", mode:'copy'
  
    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(db)

    output:
    tuple val(sample_id), path("${sample_id}*.bam"), emit: alignment
    tuple val(sample_id), path("${sample_id}*.bam.bai"), emit: alignmentindex, optional: true
    tuple val(sample_id), path("${sample_id}*report.tsv"), emit: reports, optional: true
    tuple val(sample_id), path("${sample_id}*_consensus.fa"), emit: consensus_seqs, optional: true
    tuple val(sample_id), path("${sample_id}*consensus_seqs.fa"), emit: consensus_main
    tuple val(sample_id), path("${sample_id}*_HPAI.tsv"), emit: HPAI, optional: true
    tuple val(sample_id), path("${sample_id}*_cov.png"), emit: coverage_plot, optional: true
    tuple val(sample_id), path("${sample_id}*_variants.vcf"), emit: vcf, optional: true
    tuple val(sample_id), path("logs"), emit: fluviewer_logs
    tuple val(sample_id), path("${sample_id}_fluviewer_provenance.yml"), emit: provenance
    tuple val(sample_id), path("${sample_id}*_mapping_refs.fa"), emit: ref_seqs_for_mapping, optional: true
    tuple val(sample_id), path("${sample_id}_contigs_blast.tsv"), emit: contig_blast_results, optional: true
    //tuple val(sample_id), path("${sample_id}_fluviewer/spades_output"), emit: spades_results, optional: true
    tuple val(sample_id), path("${sample_id}*.png"), emit: depth_cov_plot, optional: true

    script:
    garbage_collection = params.keep_interfiles ? '-g' : ''
    OUTPATH="${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}"

    """
    printf -- "- process_name: fluviewer\\n"                   >> ${sample_id}_fluviewer_provenance.yml
    printf -- "  tools:\\n"                                    >> ${sample_id}_fluviewer_provenance.yml
    printf -- "    - tool_name: fluviewer\\n"                  >> ${sample_id}_fluviewer_provenance.yml
    printf -- "      tool_version: \$(fluviewer --version)\\n" >> ${sample_id}_fluviewer_provenance.yml
    printf -- "  databases:\\n"                                >> ${sample_id}_fluviewer_provenance.yml
    printf -- "    - database_name: ${db}\\n"                  >> ${sample_id}_fluviewer_provenance.yml
    printf -- "      database_path: \$(readlink -f ${db})\\n"  >> ${sample_id}_fluviewer_provenance.yml
    printf -- "      database_sha256: \$(shasum -a 256 ${db}|awk '{print \$1}')\\n" >> ${sample_id}_fluviewer_provenance.yml
  
    EXITCODE=0
    (fluviewer \
	--threads ${task.cpus} \
	--forward-reads ${reads_1} \
	--reverse-reads ${reads_2} \
	--outdir . \
	--output-name ${sample_id} \
	--db ${db} \
	--min-depth ${params.min_depth}  \
	--min-mapping-quality ${params.min_q} \
	--min-identity ${params.min_ident} \
	--max-memory 40 \
	--disable-garbage-collection \
	--force && EXITCODE=\$?) \
	|| EXITCODE=\$?

    function SAFE_EXIT {
        EXITCODE=\$1
        OUTPATH=\$2

        echo \$EXITCODE > .exitcode
        cp .command.* \$OUTPATH
        cp .exitcode \$OUTPATH
        exit \$EXITCODE
    }

    OUTPATH=${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}/${sample_id}

    if [[ \$OUTPATH != /* ]]; then              # catch case where params.outdir is relative path and fix OUTPATH variable
        OUTPATH=${workflow.launchDir}/\$OUTPATH
    fi


    if [ \$EXITCODE -ne 0 ]; then 
        echo "fluviewer exited with non-zero exit code. Skipping remaining analyses."
        SAFE_EXIT \$EXITCODE \$OUTPATH
    fi

    echo "Extracting NA and HA consensus sequences..."


    if [ `grep "|HA|" ${sample_id}*consensus_seqs.fa` ]; then 
        grep -A1 "|HA|" ${sample_id}*consensus_seqs.fa > ${sample_id}_HA_consensus.fa
    else
        echo "No HA consensus sequence generated."
    fi

    if [ `grep "|NA|" ${sample_id}*consensus_seqs.fa` ]; then 
        grep -A1 "|NA|" ${sample_id}*consensus_seqs.fa > ${sample_id}_NA_consensus.fa
    else
        echo "No NA consensus sequence generated."
    fi

    if [[ ! -f ${sample_id}_HA_consensus.fa ]]; then
        echo "HA segment consensus not generated. Skipping FindCleave.py..."
    else
	python ${projectDir}/bin/FindCleave.py -i ${sample_id}_HA_consensus.fa -o ${sample_id}_HPAI.tsv
        echo "Finished running FindCleave.py."
    fi

    cp analysis_by_stage/02_blast_contigs/${sample_id}_contigs_blast.tsv .

    SAFE_EXIT \$EXITCODE \$OUTPATH
    """
}