 #!/usr/bin/env nextflow

/*
 *   A nextflow wrapper for running FluViewer
 *   -----------------------------------------
 */

nextflow.enable.dsl = 2

include { hash_files }          from './modules/hash_files.nf'
include { pipeline_provenance } from './modules/provenance.nf'
include { collect_provenance }  from './modules/provenance.nf'
include { fastp }               from './modules/fastp.nf'
include { cutadapt}             from './modules/cutadapt.nf'
include { sample_reads }     	from './modules/assembly.nf'
include { assembly }           	from './modules/assembly.nf'
include { multiqc }             from './modules/multiqc.nf'
include { fastqc }              from './modules/fastqc.nf'
include { clade_calling }       from './modules/clade_calling.nf'
include { snp_calling }         from './modules/snp_calling.nf'
include { pull_genoflu }        from './modules/genoflu.nf'
include { checkout_genoflu }    from './modules/genoflu.nf'
include { genoflu }             from './modules/genoflu.nf'

include { ref_mapping }    from './workflows/ref_mapping.nf'



// prints to the screen and to the log
log.info """
  Flu Scanner Pipeline
  ===================================
  projectDir        : ${projectDir}
  launchDir         : ${launchDir}
  database          : ${params.db}
  primers           : ${params.primers}
  fastqInputDir     : ${params.fastq_input}
  outdir            : ${params.outdir}
  pipeline run name : ${workflow.runName}
  pipeline version  : ${workflow.manifest.version}
  run_name          : ${params.run_name}
  user              : ${workflow.userName}
  Git repository    : ${workflow.repository}
  git commit id     : ${workflow.commitId}
  branch            : ${workflow.revision}
""".stripIndent()


workflow {

    ch_workflow_metadata = Channel.value([
	workflow.sessionId,
	workflow.runName,
	workflow.manifest.name,
	workflow.manifest.version,
	workflow.start,
    ])
    
    ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)

    ch_primers = Channel.fromPath(params.primer_path)

    if (params.samplesheet_input != 'NO_FILE') {
	ch_fastq_input = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }
    } else {
	ch_fastq_input = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
    }

    ch_aa_database = Channel.of([file(params.aa_database).parent, file(params.aa_database).name]).first()
    ch_nt_database = Channel.of([file(params.nt_database).parent, file(params.nt_database).name]).first()

    main:

    // Provenance channel starts with just the sample IDs
    // These will be joined to various provenance files as they are generated
    ch_provenance = ch_fastq_input.map{ it -> it[0] }

    // Generate hashes for input files
    hash_files(ch_fastq_input.map{ it -> [it[0], [it[1], it[2]]] }.combine(Channel.of("fastq_input")))

    // Clean up reads - remove adapters (fastp) and primers (cutadapt)
    fastp(ch_fastq_input)
    cutadapt(fastp.out.trimmed_reads.combine(ch_primers))
    fastqc(cutadapt.out.primer_trimmed_reads)

    sample_reads(cutadapt.out.primer_trimmed_reads)

    // Run assembly 
    assembly(sample_reads.out.reads)

    //Collect al the relevant filesfor multiqc
    ch_fastqc_collected = fastqc.out.zip.map{ it -> [it[1], it[2]]}.collect()
    multiqc(fastp.out.json.mix( cutadapt.out.log, ch_fastqc_collected ).collect().ifEmpty([]) )
 
    ref_mapping(
        ch_fastq_input,
        assembly.out.contigs,
        ch_nt_database
    )

    //Call clades for H1 and H3 samples
    //clade_calling(ref_mapping.out.ha_consensus_seq)
     
    snp_calling(ref_mapping.out.consensus, ch_aa_database)
   
    pull_genoflu(params.genoflu_github_url)

    checkout_genoflu(pull_genoflu.out.repo, params.genoflu_version)

    genoflu(ref_mapping.out.consensus.combine(pull_genoflu.out.repo))

    //
    // Provenance collection processes
    // The basic idea is to build up a channel with the following structure:
    // [sample_id, [provenance_file_1.yml, provenance_file_2.yml, provenance_file_3.yml...]]
    // ...and then concatenate them all together in the 'collect_provenance' process.
    ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it ->      [it[0], [it[1]]] }
    ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it ->      [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(fastp.out.provenance).map{ it ->           [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(cutadapt.out.provenance).map{ it ->        [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(sample_reads.out.provenance).map{ it ->    [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(assembly.out.provenance).map{ it ->        [it[0], it[1] << it[2]] }    
    //ch_provenance = ch_provenance.join(clade_calling.out.provenance).map{ it ->   [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(snp_calling.out.provenance).map{ it ->     [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(genoflu.out.provenance).map{ it ->         [it[0], it[1] << it[2]] }
    collect_provenance(ch_provenance)

}
