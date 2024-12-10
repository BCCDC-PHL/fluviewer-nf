process blastn {

    errorStrategy 'ignore'

    label 'medium'

    tag {sample_id}

    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}*blastn.tsv" , mode:'copy'
    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}*filter.tsv" , mode:'copy'
    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}*fasta" , mode:'copy'

    input: 
    tuple val(sample_id), path(contigs)        // 1) Name  2) Assembled contigs  3) FASTA file where the final output sequences are stored
    tuple path(blast_db_path), val(blast_db_name)                                                          // BLAST database used to search for best matching reference 

    output:
    tuple val(sample_id), path("${sample_id}*blastn.tsv"), emit: raw
    tuple val(sample_id), path("${sample_id}*filter.tsv"), path("${sample_id}*fasta"), emit: main
    tuple val(sample_id), path("${sample_id}*_ref.fasta"), emit: ref
    path("${sample_id}*filter.tsv"), emit: filter
    tuple val(sample_id), path("${sample_id}-*-provenance.yml"), emit: provenance

    script:
    
    """
    printf -- "- process_name: blastn\\n" > ${sample_id}-blastn-provenance.yml
    printf -- "  tools: \\n  - tool_name: blastn\\n    tool_version: \$(blastn -version 2>&1 | head -n1 | cut -d' ' -f2)\\n" >> ${sample_id}-blastn-provenance.yml

    blastn -num_threads ${task.cpus} -query ${contigs} -db ${blast_db_path}/${blast_db_name} -outfmt 6 > ${sample_id}_blastn.tsv &&

    filter_alignments.py ${sample_id}_blastn.tsv \
    --seqs ${blast_db_path}/${blast_db_name} \
    --out_tsv ${sample_id}_blastn_filter.tsv \
    --out_fasta ${sample_id}_ref.fasta 
    """
}


process map_reads {

    tag {sample_id}

    conda ""

    errorStrategy 'ignore'

	label 'heavy'

    input: 
    tuple val(sample_id), path(reads_1), path(reads_2), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.sam"), emit: sam
    tuple val(sample_id), path("${sample_id}-*-provenance.yml"), emit:provenance

    script:
	"""
    printf -- "- process_name: map_reads\\n" > ${sample_id}-bwa-provenance.yml
    printf -- "  tools: \\n  - tool_name: bwa\\n    tool_version: \$(bwa 2>&1 | head -n3 | tail -n1 | cut -d' ' -f2)\\n" >> ${sample_id}-bwa-provenance.yml

    bwa index ${reference}
	bwa mem -t ${task.cpus} ${reference} ${reads_1} ${reads_2} > ${sample_id}.sam
	"""
}
process sort_filter_sam {

	tag {sample_id}

    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}*{depth.tsv}", mode: 'copy'


    input: 
    tuple val(sample_id), path(samfile)

    output:
    tuple val(sample_id), path("${sample_id}.sort.filter.bam"), path("${sample_id}*bai"), emit: bam
    tuple val(sample_id), path("${sample_id}_depth.tsv"), emit: depth
    tuple val(sample_id), path("${sample_id}-*-provenance.yml"), emit: provenance

    script:
	"""
    printf -- "- process_name: sort_filter_sam\\n" > ${sample_id}-samtools-provenance.yml
    printf -- "  tools: \\n  - tool_name: samtools\\n    tool_version: \$(samtools --version 2>&1 | head -n1 | cut -d' ' -f2)\\n" >> ${sample_id}-samtools-provenance.yml

	samtools view -f 3 -F 2828 -q 30 -h ${samfile} | samtools sort -o ${sample_id}.sort.filter.bam 
    samtools index ${sample_id}.sort.filter.bam 

    samtools depth -a ${sample_id}.sort.filter.bam  > ${sample_id}_depth.tsv
	"""
}

process make_pileup {

	tag {sample_id}

    conda "${projectDir}/environments/fluviewer.yml"

    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}*tsv", mode: 'copy'

    input: 
    tuple val(sample_id), path(bam), path(bam_index)

    output:
    tuple val(sample_id), path("${sample_id}*pileup.tsv"), emit: pileup

    script:
	"""
    make_pileup.py --input ${bam} --output ${sample_id}_pileup.tsv
	"""
}



process freebayes {

    tag { sample_id }

    conda "${projectDir}/environments/fluviewer.yml"

    publishDir "${params.outdir}/${sample_id}/vcf", pattern: "${sample_id}*{gz,csi}", mode: 'copy'

    input: 
	tuple val(sample_id), path(bamfile), path(bam_index), path(reference)
	
	output:
    tuple val(sample_id), path("${sample_id}*vcf.gz"), path("${sample_id}*.csi"), emit: vcf
    tuple val(sample_id), path("${sample_id}*provenance.yml"), emit: provenance

    script:

    """
	printf -- "- process_name: freebayes\\n"                                                                                           >> ${sample_id}_freebayes_provenance.yml
    printf -- "  tools: \\n  - tool_name: samtools\\n    tool_version: \$(samtools --version 2>&1 | head -n1 | cut -d' ' -f2)\\n"   >> ${sample_id}_freebayes_provenance.yml
    printf -- "  - tool_name: freebayes\\n    tool_version: \$(freebayes --version | cut -d' ' -f3)\\n"                             >> ${sample_id}_freebayes_provenance.yml
    printf -- "  - tool_name: bcftools\\n    tool_version: \$(bcftools --version | head -n1 | cut -d' ' -f2)\\n"                    >> ${sample_id}_freebayes_provenance.yml

    samtools faidx ${reference}
    freebayes -X -O \
    -b ${bamfile} \
    -v ${sample_id}.vcf \
    -f ${reference} \
    --limit-coverage 200 \
    -p 1 \
    -m 20 \
    -q 20 \
    --pooled-continuous \

    bcftools sort --output-type z ${sample_id}.vcf > ${sample_id}.vcf.gz
    bcftools index ${sample_id}.vcf.gz
    """

}

process snpeff_annotation {
    tag { sample_id }

    conda "${projectDir}/environments/fluviewer.yml"

    publishDir "${params.outdir}/${sample_id}/vcf", pattern: "${sample_id}*{vcf}", mode: 'copy'

    input: 
	tuple val(sample_id), path(vcf), path(vcf_index), path(snpeff_db_path)
	
	output:
    tuple val(sample_id), path("${sample_id}*vcf"), emit: vcf
    tuple val(sample_id), path("${sample_id}*provenance.yml"), emit: provenance

    script:

    """
	printf -- "- process_name: snpeff_annotation\\n"                                                      >> ${sample_id}_snpeff_provenance.yml
    printf -- "  tools: \\n  "                                                                            >> ${sample_id}_snpeff_provenance.yml
    printf -- "  - tool_name: snpeff\\n    tool_version: \$(snpEff -version 2>&1 | cut -d\$'\t' -f2)\\n"   >> ${sample_id}_snpeff_provenance.yml
    printf -- "  - tool_name: samtools\\n    tool_version: \$(samtools --version 2>&1 | head -n1 | cut -d' ' -f2)\\n"   >> ${sample_id}_snpeff_provenance.yml
    printf -- "  - tool_name: bcftools\\n    tool_version: \$(bcftools --version | head -n1 | cut -d' ' -f2)\\n"   >> ${sample_id}_snpeff_provenance.yml

    run_snpeff.py --input ${vcf} --config ${snpeff_db_path} --outname ${sample_id}
    """

}


process mask_low_coverage {
	tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}*.bed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bamfile), path(bam_index)
	
    output:
    tuple val(sample_id), path("${sample_id}_low_coverage.bed"), emit: masked_bed
    tuple val(sample_id), path("${sample_id}-*-provenance.yml"), emit: provenance


    script:
    """
    printf -- "- process_name: mask_low_coverage\\n" > ${sample_id}-bedtools-provenance.yml
    printf -- "  tools: \\n  - tool_name: bedtools\\n    tool_version: \$(bedtools --version 2>&1 | cut -d' ' -f2)\\n" >> ${sample_id}-bedtools-provenance.yml

    bedtools genomecov -bga -ibam ${bamfile} |  awk '\$4 < ${params.min_depth} {{print}}' | awk 'BEGIN{FS=OFS="\\t"} {print \$1,\$2+1,\$3+1,\$4}' > ${sample_id}_low_coverage.bed
    """
}

process make_consensus {

    tag { sample_id }

    conda "${projectDir}/environments/fluviewer.yml"

    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}*consensus.fasta", mode: 'copy'
    
    input:
    tuple val(sample_id), path(vcf), path(vcf_index), path(reference), path(mask_file)
	
    output:
    tuple val(sample_id), path("${sample_id}*consensus.fasta"), emit: consensus
    tuple val(sample_id), path("${sample_id}-*-provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: make_consensus\\n" > ${sample_id}-bcftools-provenance.yml
    printf -- "  tools: \\n  - tool_name: bcftools\\n    tool_version: \$(bcftools --version 2>&1 | head -n1 | cut -d' ' -f2)\\n" >> ${sample_id}-bcftools-provenance.yml

    bcftools consensus -m ${mask_file} -f ${reference} ${vcf} > ${sample_id}_consensus.fasta
    """
}