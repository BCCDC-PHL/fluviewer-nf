process normalize_depth {
    tag { sample_id }

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}-normalized_R1.fastq.gz"), path("${sample_id}-normalized_R2.fastq.gz"), emit: normalized_reads
    tuple val(sample_id), path("${sample_id}_normalize_depth_provenance.yml"), emit: provenance

    script:
    max_memory_gb = task.memory.toString().split(" ")[0]
    """
    printf -- "- process_name: normalize_depth\\n"                >> ${sample_id}_normalize_depth_provenance.yml
    printf -- "  tools:\\n"                                       >> ${sample_id}_normalize_depth_provenance.yml
    printf -- "    - tool_name: bbnorm\\n"                        >> ${sample_id}_normalize_depth_provenance.yml
    printf -- "      tool_version: \$(bbnorm.sh --version 2>&1 | head -n 2 | tail -n 1 | cut -d ' ' -f 3)\\n" >> ${sample_id}_normalize_depth_provenance.yml
    
    bbnorm.sh \
	-Xmx${max_memory_gb}g \
	in1=${reads_1} \
	in2=${reads_2} \
	out1=${sample_id}-normalized_R1.fastq \
	out2=${sample_id}-normalized_R2.fastq \
	target=${params.target_depth} \
	

    gzip ${sample_id}-normalized_R1.fastq
    gzip ${sample_id}-normalized_R2.fastq
    """
    
}

process fluviewer {

    tag { sample_id }

    errorStrategy 'ignore'

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}*", mode:'copy', saveAs: { filename -> filename.split("/").last() }
    publishDir "${params.outdir}/${sample_id}", pattern: "*tsv", mode:'copy', saveAs: { filename -> filename.split("/").last() }
    publishDir "${params.outdir}/${sample_id}", pattern: "logs", mode:'copy', saveAs: { filename -> "fluviewer_logs" }  
  
    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(db)

    output:
    tuple val(sample_id), path("${sample_id}*.bam"), emit: alignment
    tuple val(sample_id), path("${sample_id}*.bam.bai"), emit: alignmentindex, optional: true
    tuple val(sample_id), path("${sample_id}*report.tsv"), emit: reports, optional: true
    tuple val(sample_id), path("${sample_id}*_consensus.fa"), emit: consensus_seqs, optional: true
    tuple val(sample_id), path("${sample_id}_HA_consensus.fa"), emit: ha_consensus_seq, optional: true
    tuple val(sample_id), path("${sample_id}*consensus_seqs.fa"), emit: consensus_main
    tuple val(sample_id), path("${sample_id}*_HPAI.tsv"), emit: HPAI, optional: true
    tuple val(sample_id), path("${sample_id}*_cov.png"), emit: coverage_plot, optional: true
    tuple val(sample_id), path("${sample_id}*_variants.vcf"), emit: vcf, optional: true
    tuple val(sample_id), path("logs"), emit: fluviewer_logs
    tuple val(sample_id), path("${sample_id}_fluviewer_provenance.yml"), emit: provenance
    tuple val(sample_id), path("${sample_id}*_mapping_refs.fa"), emit: ref_seqs_for_mapping, optional: true
    tuple val(sample_id), path("${sample_id}_contigs_blast.tsv"), emit: contig_blast_results, optional: true
    tuple val(sample_id), path("${sample_id}*.png"), emit: depth_cov_plot, optional: true

    script:
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
	--skip-depth-normalization \
	--force && EXITCODE=\$?) \
	|| EXITCODE=\$?

    function SAVE_LOGS {
        EXITCODE=\$1
        
        mkdir -p logs
        echo \$EXITCODE > logs/fluviewer_exitcode.txt
        cp .command.out logs/fluviewer_stdout.txt
        cp .command.err logs/fluviewer_stderr.txt
    }

    function SAFE_EXIT {
        EXITCODE=\$1
        OUTPATH=\$2

        mkdir -p \${OUTPATH}/fluviewer_logs
        cp logs/fluviewer*.txt \${OUTPATH}/fluviewer_logs
        exit \$EXITCODE
    }

    OUTPATH=${params.outdir}/${sample_id}

    if [[ \$OUTPATH != /* ]]; then              # catch case where params.outdir is relative path and fix OUTPATH variable
        OUTPATH=${workflow.launchDir}/\$OUTPATH
    fi

    SAVE_LOGS \$EXITCODE

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
        FindCleave.py -i ${sample_id}_HA_consensus.fa -o ${sample_id}_HPAI.tsv
        echo "Finished running FindCleave.py."
    fi

    cp analysis_by_stage/*_blast_contigs/${sample_id}_contigs_blast.tsv .

    """
}
