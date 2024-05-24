 #!/usr/bin/env nextflow

/*
 *   A nextflow wrapper for running FluViewer
 *   -----------------------------------------

 == V1  ==
This pipeline will run FluViewer on a set of fastq files in a baseDir.
Each output will be its own directory.
Future versions will add in:
- fastp to remove adapters and produce QC
- multqc to read the results of this
- a script to scrape together results and produce an output csv

 */

 import java.time.LocalDateTime

 nextflow.enable.dsl = 2

 include { hash_files } from './modules/hash_files.nf'
 include { pipeline_provenance } from './modules/provenance.nf'
 include { collect_provenance } from './modules/provenance.nf'
 include { fastp } from './modules/fastp.nf'
 include { cutadapt} from './modules/cutadapt.nf'
 include { FluViewer } from './modules/FluViewer.nf'
 include { multiqc } from './modules/multiqc.nf'
 include { FASTQC } from './modules/fastqc.nf'
 include { CLADE_CALLING } from './modules/clade_calling.nf'
 include { SNP_CALLING } from './modules/snp_calling.nf'
 include { PULL_GENOFLU ; CHECKOUT_GENOFLU ;  GENOFLU } from './modules/genoflu.nf'


// prints to the screen and to the log
        log.info """

                 FluViewer Pipeline
                 ===================================
                 projectDir        : ${projectDir}
                 launchDir         : ${launchDir}
                 database          : ${params.db}
                 primers           : ${params.primers}
                 fastqInputDir     : ${params.fastq_input}
                 outdir            : ${params.outdir}
                 pipeline run      : ${params.pipeline_short_name}
                 pipeline version  : ${params.pipeline_minor_version}
                 run_name          : ${params.run_name}
                 user              : $workflow.userName
                 Git repository    : $workflow.repository
                 git commit id     : $workflow.commitId
                 branch            : $workflow.revision
                 """
                 .stripIndent()


workflow {
     ch_start_time = Channel.of(LocalDateTime.now())
     ch_pipeline_name = Channel.of(workflow.manifest.name)
     ch_pipeline_version = Channel.of(workflow.manifest.version)


     ch_pipeline_provenance = pipeline_provenance(ch_pipeline_name.combine(ch_pipeline_version).combine(ch_start_time))

     ch_primers = Channel.fromPath(params.primer_path)
     ch_db = Channel.fromPath(params.db)
     ch_fastq_input = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }

     ch_reference_db = Channel.of([file(params.blastx_subtype_db).parent, file(params.blastx_subtype_db).name]).first()

     main:
     // Generate hashes for input files
     hash_files(ch_fastq_input.map{ it -> [it[0], [it[1], it[2]]] }.combine(Channel.of("fastq_input")))

     // Clean up reads - remove adapters (fastp) and primers (cutadapt)
     fastp( ch_fastq_input )
     cutadapt(fastp.out.trimmed_reads.combine(ch_primers))
     FASTQC(cutadapt.out.primer_trimmed_reads)
     
     // Run FluViewer 
     FluViewer(cutadapt.out.primer_trimmed_reads.combine(ch_db))

     //Collect al the relevant filesfor MULTIQC
     ch_fastqc_collected = FASTQC.out.zip.map{ it -> [it[1], it[2]]}.collect()
     multiqc(fastp.out.json.mix( cutadapt.out.log, ch_fastqc_collected ).collect().ifEmpty([]) )
 
     //Call clades for H1 and H3 samples
     CLADE_CALLING(FluViewer.out.consensus_seqs)
     
     SNP_CALLING(FluViewer.out.consensus_main, ch_reference_db)

     
     PULL_GENOFLU(params.genoflu_github_url)

     CHECKOUT_GENOFLU(PULL_GENOFLU.out.repo, params.genoflu_version)

     GENOFLU(FluViewer.out.consensus_main.combine(PULL_GENOFLU.out.repo))



     //Pool Provenance data
     ch_provenance = FluViewer.out.provenance
     ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it -> [it[0], [it[1]] << it[2]] }
     ch_provenance = ch_provenance.join(fastp.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
     ch_provenance = ch_provenance.join(cutadapt.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
     ch_provenance = ch_provenance.join(ch_fastq_input.map{ it -> it[0] }.combine(ch_pipeline_provenance)).map{ it -> [it[0], it[1] << it[2]] }
     ch_provenance = ch_provenance.join(CLADE_CALLING.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
     ch_provenance = ch_provenance.join(GENOFLU.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
     collect_provenance(ch_provenance)



}