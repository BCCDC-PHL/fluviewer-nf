 #!/usr/bin/env nextflow

/*
 *   A nextflow wrapper for running FluViewer
 *   -----------------------------------------
 */

nextflow.enable.dsl = 2


include { map_reads }           from '../modules/ref_mapping.nf'
include { sort_filter_sam }     from '../modules/ref_mapping.nf'
include { make_pileup }         from '../modules/ref_mapping.nf'
include { blastn }              from '../modules/ref_mapping.nf'
include { freebayes }           from '../modules/ref_mapping.nf'
include { make_coverage_plot }  from '../modules/plotting.nf'
include { mask_low_coverage }   from '../modules/ref_mapping.nf'
include { make_consensus }      from '../modules/ref_mapping.nf'
include { snpeff_annotation }   from '../modules/ref_mapping.nf'
include { mutation_watch    }   from '../modules/ref_mapping.nf'



workflow ref_mapping {

    take:
        ch_fastq
        ch_contigs
        ch_nt_database

    main:
        blastn(ch_contigs, ch_nt_database)
        map_reads(ch_fastq.join(blastn.out.ref))
        sort_filter_sam(map_reads.out.sam)
        make_pileup(sort_filter_sam.out.bam)
        make_coverage_plot(sort_filter_sam.out.depth)
        freebayes(sort_filter_sam.out.bam.join(blastn.out.ref))

        if (params.snpeff_config != 'NO_FILE') {
            snpeff_annotation(freebayes.out.vcf.combine(Channel.fromPath(params.snpeff_config)))
            mutation_watch(snpeff_annotation.out.vcf_filter.combine(Channel.fromPath(params.mutation_watchlist)))
        }

        mask_low_coverage(sort_filter_sam.out.bam)
        make_consensus(freebayes.out.vcf.join(blastn.out.ref).join(mask_low_coverage.out.masked_bed))
    emit:
        consensus = make_consensus.out.consensus
}
