process cutadapt {

    tag { sample_id }

    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(primers)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trim_p.fastq.gz"), path("${sample_id}_R2.trim_p.fastq.gz"), emit: primer_trimmed_reads
    path("${sample_id}.cutadapt.log"), emit: log
    tuple val(sample_id), path("${sample_id}_cutadapt_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: cutadapt\\n"    >> ${sample_id}_cutadapt_provenance.yml
    printf -- "  tools:\\n"                    >> ${sample_id}_cutadapt_provenance.yml
    printf -- "    - tool_name: cutadapt\\n"   >> ${sample_id}_cutadapt_provenance.yml
    printf -- "      tool_version: \$(cutadapt --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_cutadapt_provenance.yml

    cutadapt \
      -j ${task.cpus} \
      -a file:${params.primers} \
      -A file:${params.primers} \
      -g file:${params.primers} \
      -G file:${params.primers} \
      -o ${sample_id}_R1.trim_p.fastq.gz \
      -p ${sample_id}_R2.trim_p.fastq.gz \
      ${sample_id}_R1.trim.fastq.gz \
      ${sample_id}_R2.trim.fastq.gz \
      > ${sample_id}.cutadapt.log
    """
}
