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
check_param_list = []
for (param in check_param_list) { if (param) { file(param, checkIfExists: true) } }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DORADO_BASECALLER } from './modules/local/dorado/basecaller/main'
include { NANOPLOT } from './modules/nf-core/nanoplot/main'
include { PYCOQC } from './modules/nf-core/pycoqc/main'
include { MULTIQC } from './modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    //
    // INIT: 
    // 
    ch_versions = Channel.empty{}


    //
    // CHANNEL: Adding all pod5 files
    //
    ch_pod5_files_pass = Channel.fromPath("${params.run_dir}/pod5_pass/*.pod5")
    ch_pod5_files_fail = Channel.fromPath("${params.run_dir}/pod5_fail/*.pod5")
    ch_pod5_files_skipped = Channel.fromPath("${params.run_dir}/pod5_skipped/*.pod5")
    ch_pod5_files = ch_pod5_files_pass.mix(ch_pod5_files_fail).mix(ch_pod5_files_skipped)

    

    // TODO: use run ID as sample name
    // ch_pod5_files = ch_pod5_files
    //     .collect()
    //     .map{
    //         [[ id: it.baseName ], [ it ] ]
    //     }

    ch_pod5_files = ch_pod5_files
        .collect()
        .map{
            [[ id: 'RUN' ], it ]
        }


    // ch_pod5_files | view

    // Include fastq data and link it to its metadata in a channel
    // ch_fastq = Channel.fromPath("${params.fastq}/*.fastq*")
    // ch_txt = Channel.fromPath("${params.txt}/*.txt")

    //
    // CHANNEL: Adding metadata
    //

    // ch_fastq = ch_fastq.map {
    //     [ [id: it.baseName ], [ it ] ]
    // }

    // ch_txt = ch_txt.map {
    //     [ [id: it.baseName ], [ it ] ]
    // }

    // ch_txt | view

    //
    // MODULE: Run nanoplot, pycoqc
    // 

    DORADO_BASECALLER (
        ch_pod5_files,
        params.dorado_model
    )

    // NANOPLOT (
    //     ch_fastq
    // )

    // PYCOQC (
    //     ch_txt
    // )


    // 
    // MODULE: MULTIQC
    // 

    // workflow_summary = multiqc_summary(workflow, params)
    // ch_workflow_summary = Channel.value(workflow_summary)

    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // // ch_multiqc_files = ch_multiqc_files.mix(DUMP_SOFTWARE_VERSIONS.out.mqc_yml.collect())
    // // ch_multiqc_files = ch_multiqc_files.mix(DUMP_SOFTWARE_VERSIONS.out.mqc_unique_yml.collect())

    // // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    // // ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT_FASTQ.out.txt.collect{it[1]}.ifEmpty([]))

    // MULTIQC (
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config,
    //     [],
    // //     []
    // )
    // multiqc_report = MULTIQC.out.report.toList()

    // ch_versions  = ch_versions.mix(NANOPLOT.out.versions)

}
