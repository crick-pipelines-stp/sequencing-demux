#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pipelines-technology/nanopore_demux
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/crick-pipelines-stp/nanopore-demux
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { summary_log     } from './modules/local/util/logging/main'
include { multiqc_summary } from './modules/local/util/logging/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INIT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

log.info summary_log(workflow, params, params.debug, params.monochrome_logs)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check manditory input parameters to see if the files exist if they have been specified
check_param_list = [
    run_dir: params.run_dir
]
for (param in check_param_list) {
    if (!param.value) {
        exit 1, "Required parameter not specified: ${param.key}"
    }
    else {
        file(param.value, checkIfExists: true)
    }
}

// Check non-manditory input parameters to see if the files exist if they have been specified
check_param_list = [
    params.bam,
    params.samplesheet
]
for (param in check_param_list) { if (param) { file(param, checkIfExists: true) } }

// Select dorado model
// Dorado models can be selected by auto selection (e.g. hac selects the latest compatible hac model) or by direct model selection.
// In the case of direct selection, we need to add the container path to the model.
dorado_model = params.dorado_auto_model
if(params.dorado_model) {
    dorado_model = "/home/" + params.dorado_model
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DORADO_BASECALLER        } from './modules/local/dorado/basecaller/main'
include { DORADO_DEMUX             } from './modules/local/dorado/demux/main'
include { SAMTOOLS_VIEW            } from './modules/nf-core/samtools/view/main'
include { CHOPPER                  } from './modules/nf-core/chopper/main'
include { FASTQC                   } from './modules/nf-core/fastqc/main'
include { NANOPLOT                 } from './modules/nf-core/nanoplot/main'
include { PYCOQC                   } from './modules/nf-core/pycoqc/main'
include { MULTIQC                  } from './modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_USER  } from './modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMPLESHEET_PARSE } from './subworkflows/local/samplesheet_parse'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    //
    // INIT:
    // 
    ch_versions    = Channel.empty()
    ch_samplesheet = Channel.empty()

    //
    // CHANNEL: Adding all pod5 files
    //
    ch_pod5_files_pass    = Channel.fromPath("${params.run_dir}/pod5_pass/*.pod5")
    ch_pod5_files_fail    = Channel.fromPath("${params.run_dir}/pod5_fail/*.pod5")
    ch_pod5_files_skipped = Channel.fromPath("${params.run_dir}/pod5_skipped/*.pod5")
    ch_pod5_files = ch_pod5_files_pass.mix(ch_pod5_files_fail).mix(ch_pod5_files_skipped)

    //
    // CHANNEL: Put all pod5 generated files and their corresponding sample IDs into a single channel 
    //
    ch_pod5_files = ch_pod5_files
        .collect()
        .map{ [[ id: 'pod5' ], it ] }

    //
    // CHANNEL: Search and add sequencing summary
    //
    ch_sequencing_summary = Channel.fromPath("${params.run_dir}/sequencing_summary*.txt", checkIfExists: true)
        .map{ [ [ id: it.simpleName ], it ] }

    //
    // CHANNEL: Adding bam files to a channel if it exists
    //
    ch_bam = Channel.empty()
    if (params.bam) {
        ch_bam = Channel.from(file(params.bam, checkIfExists: true))
            .map{ [ [ id: it.simpleName ], it ] }
    }

    //
    // CHANNEL: Load samplesheet
    //
    if(params.samplesheet) {
        ch_samplesheet = Channel.from(file(params.samplesheet))
    }

    // 
    // SUBWORKFLOW: check input samplesheet and add relevant info to metadata
    // 
    SAMPLESHEET_PARSE (
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(SAMPLESHEET_PARSE.out.versions)
    ch_meta     = SAMPLESHEET_PARSE.out.meta

    //
    // CHANNEL: extract run ID name and assign to metadata
    //
    runid = file(params.run_dir).name
    ch_meta = ch_meta.map{
        it.run_id = runid 
        it.id = it.sample_id
        it.remove("sample_id")
        it 
    }
    
    //
    // MODULE: Generate a bam file using the Dorado basecaller unless a bam file was already present as an input
    // 
    if (!params.bam) {
        DORADO_BASECALLER (
            ch_pod5_files,
            dorado_model
        )
        ch_versions = ch_versions.mix(DORADO_BASECALLER.out.versions)
        ch_bam      = DORADO_BASECALLER.out.bam
    }

    // 
    // MODULE: Generate demultiplexed bam or fastq files
    // 
    DORADO_DEMUX (
        ch_bam
    )
    ch_versions    = ch_versions.mix(DORADO_DEMUX.out.versions)
    ch_demux_bam   = DORADO_DEMUX.out.bam
    ch_demux_fastq = DORADO_DEMUX.out.fastq

    //
    // CHANNEL: Merge metadata to the demultiplexed file
    //
    ch_demux_fastq = ch_meta
        .map { [it.barcode, it] }
        .join( ch_demux_fastq.map{ [ it[1].simpleName, it ] } )
        .map { [ it[1], it[2][1] ] }

    // ch_demux_fastq |view

    //
    // CHANNEL: Merge metadata to the demultiplexed file
    //
    ch_demux_bam = ch_meta
        .map { [it.barcode, it] }
        .join( ch_demux_bam.map{ [ it[1].simpleName, it ] } )
        .map { [ it[1], it[2][1] ] }

    // 
    // MODULE: Run chopper on fastq files
    //
    CHOPPER (
        ch_demux_fastq
    )
    ch_versions      = ch_versions.mix(CHOPPER.out.versions)
    ch_chopper_fastq = CHOPPER.out.fastq

    // 
    // MODULE: Run fastqc on fastq files
    //
    FASTQC (
        ch_chopper_fastq
    )
    ch_versions      = ch_versions.mix(FASTQC.out.versions)
    ch_fastqc_zip = FASTQC.out.zip

    // 
    // MODULE: Filter reads
    // 
    SAMTOOLS_VIEW (
        ch_demux_bam,
        [],
        []
    )
    ch_versions    = ch_versions.mix(SAMTOOLS_VIEW.out.versions)
    ch_filtered_reads_bam   = SAMTOOLS_VIEW.out.bam

    //
    // MODULE: Run Nanoplot
    //
    // NANOPLOT (
    //     ch_sequencing_summary
    // )
    // ch_versions = ch_versions.mix(NANOPLOT.out.versions)
    
    // 
    // MODULE: Run PYCOQC
    //
    PYCOQC (
        ch_sequencing_summary
    )
    ch_versions = ch_versions.mix(PYCOQC.out.versions)
    ch_pycoqc   = PYCOQC.out.json

    // 
    // MODULE: MULTIQC
    // 
    workflow_summary = multiqc_summary(workflow, params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(DUMP_SOFTWARE_VERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(DUMP_SOFTWARE_VERSIONS.out.mqc_unique_yml.collect())

    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT.out.txt.collect{it[1]}.ifEmpty([]))

    ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_pycoqc.collect{it[1]}.ifEmpty([]))

    // ch_multiqc_files |view
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config,
        [],
        []
    )
    multiqc_report = MULTIQC.out.report.toList()


    // ch_multiqc_user_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_user_files = ch_multiqc_user_files.mix(ch_fastqc_zip.collect{it[1]}.ifEmpty([]))
    // // ch_multiqc_user_files = ch_multiqc_user_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    // // ch_multiqc_user_files | view

    // MULTIQC_USER (
    //     ch_multiqc_user_files.collect(),
    //     ch_multiqc_config,
    //     [],
    //     []
    // )
    // multiqc_user_report = MULTIQC_USER.out.report.toList()
}
