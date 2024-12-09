process sample_reads {
    tag { sample_id }

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_R1_sample.fq.gz"), path("${sample_id}_R2_sample.fq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}_sample_reads_provenance.yml"), emit: provenance

    script:
    max_memory_gb = task.memory.toString().split(" ")[0]
    """
    printf -- "- process_name: sample_reads\\n"                >> ${sample_id}_sample_reads_provenance.yml
    printf -- "  tools:\\n"                                       >> ${sample_id}_sample_reads_provenance.yml
    printf -- "    - tool_name: seqtk\\n"                        >> ${sample_id}_sample_reads_provenance.yml
    printf -- "      tool_version: \$(seqtk 2>&1 | grep "Version" | cut -d' ' -f2)\\n" >> ${sample_id}_sample_reads_provenance.yml
    
    READ_COUNT_1=\$(zcat ${reads_1} | grep "^@"  | wc -l)
    READ_COUNT_2=\$(zcat ${reads_2} | grep "^@"  | wc -l)
    let "TARGET_COUNT = ${params.target_depth} * ${params.genome_size}"

    if [ \${READ_COUNT_1} -gt \${TARGET_COUNT} ]; then 
        seqtk sample -s 12 ${reads_1} \${TARGET_COUNT} > ${sample_id}_R1_sample.fq
        seqtk sample -s 12 ${reads_2} \${TARGET_COUNT} > ${sample_id}_R2_sample.fq

        gzip ${sample_id}_R1_sample.fq
        gzip ${sample_id}_R2_sample.fq
        echo "Completed random downsample on reads.\nOriginal read count: \${READ_COUNT_1} \nNew read count: \${TARGET_COUNT}"
    else
        ln -s ${reads_1} ${sample_id}_R1_sample.fq.gz
        ln -s ${reads_2} ${sample_id}_R2_sample.fq.gz
        echo "WARNING: Skipping random downsample step. Number of reads in original file exceed target read count. \nOriginal read count: \${READ_COUNT_1} \nAttempted target count: \${TARGET_COUNT}"
    fi

    """
    
}

process assembly {
    
    tag { sample_id } 

    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}*fasta", mode: 'copy'
    
    input:
    tuple val(sample_id), path(fwd_reads), path(rev_reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_spades"), emit: spades_dir
    tuple val(sample_id), path("${sample_id}_contigs.fasta"), emit: contigs
    tuple val(sample_id), path("${sample_id}*.yml"), emit: provenance
    
    script:
    """
    printf -- "- process_name: assembly\\n"                >> ${sample_id}_assembly_provenance.yml
    printf -- "  tools:\\n"                                       >> ${sample_id}_assembly_provenance.yml
    printf -- "    - tool_name: assembly\\n"                        >> ${sample_id}_assembly_provenance.yml
    printf -- "      tool_version: \$(spades.py --version | cut -d' ' -f4)\\n" >> ${sample_id}_assembly_provenance.yml
    
    spades.py \
        --threads ${task.cpus} \
        --rnaviral \
        --isolate \
        -1 ${fwd_reads} \
        -2 ${rev_reads} \
        -o ${sample_id}_spades
    
    cp ${sample_id}_spades/contigs.fasta ./${sample_id}_contigs.fasta
    """
}