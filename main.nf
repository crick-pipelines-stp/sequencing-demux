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

include { params_summary_map        } from './modules/local/util/logging/main'
include { summary_log               } from './modules/local/util/logging/main'
include { multiqc_summary           } from './modules/local/util/logging/main'
include { dump_parameters           } from './modules/local/util/logging/main'
include { im_notification           } from './modules/local/util/logging/main'
include { dump_meta                 } from './modules/local/util/logging/main'
include { gen_id                    } from './modules/local/util/logging/main'
include { workflow_complete_summary } from './modules/local/util/logging/main'

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
def summary_params = params_summary_map(workflow, params, params.debug)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check manditory input parameters to see if the files exist if they have been specified
check_param_list = [
    run_dir: params.run_dir,
    sample_sheet: params.samplesheet
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
    params.bam
]
for (param in check_param_list) { if (param) { file(param, checkIfExists: true) } }

// Select dorado model
// Dorado models can be selected by auto selection (e.g. hac selects the latest compatible hac model) or by direct model selection.
// In the case of direct selection, we need to add the container path to the model.
dorado_model = params.dorado_auto_model
if(params.dorado_model) {
    dorado_model = "/home/" + params.dorado_model
}

// Check bc-kit
if(params.dorado_bc_kit && !(params.dorado_bc_kit in params.bc_kits)) {
    exit 1, "Invalid barcode kit specified: ${params.dorado_bc_kit}"
}
// Extract barcode kit value from the summary file
if !(params.dorado_bc_kit) {
    // Find the summary file
    def summaryFileDir = file("${projectDir}")
    def summaryFile = summaryFileDir.list().find { it.contains('sequencing_summary') && it.endsWith('.txt') }
    // Check if a summary file was found
    if (!summaryFile) {
        exit 1, "No summary file found in ${summaryFileDir}."
    }
    // Read the entire content of the summary file as a single string
    def summaryContent = summaryFile.text

    // Find the first matching barcode kit in the summary file
    def extrapolatedBcKit = params.bc_kits.find { bc_kit ->
        summaryContent.contains(bc_kit)
    }

    // Set `params.dorado_bc_kit` to the found kit or null if no match
    params.dorado_bc_kit = extrapolatedBcKit ?: null
    if (!params.dorado_bc_kit) {
        exit 1, "No valid barcode kit found in the summary file."
    } 
}

// Extract run_id
def runid = file(params.run_dir).name

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DORADO_BASECALLER                   } from './modules/local/dorado/basecaller/main'
include { SAMTOOLS_MERGE as MERGE_BASECALLING } from './modules/nf-core/samtools/merge/main'
include { DORADO_DEMUX                        } from './modules/local/dorado/demux/main'
include { SAMTOOLS_VIEW as FILTER_READ_Q      } from './modules/nf-core/samtools/view/main'
include { CHOPPER                             } from './modules/nf-core/chopper/main'
include { FASTQC                              } from './modules/nf-core/fastqc/main'
include { SAMTOOLS_VIEW as FILTER_BC          } from './modules/nf-core/samtools/view/main'
include { SAMTOOLS_MERGE as MERGE_GROUPS      } from './modules/nf-core/samtools/merge/main'
include { DORADO_SUMMARY                      } from './modules/local/dorado/summary/main'
include { TOULLIGQC as TOULLIGQC_ALL          } from './modules/local/toulligqc/main'
include { TOULLIGQC as TOULLIGQC_GROUPED      } from './modules/local/toulligqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS         } from './modules/local/custom_dumpsoftwareversions'
include { MULTIQC as MULTIQC_ALL              } from './modules/local/multiqc/main'
include { MULTIQC as MULTIQC_GROUPED          } from './modules/local/multiqc/main'

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
    ch_pod5_files         = Channel.fromPath("${params.run_dir}/pod5/*.pod5")
    ch_pod5_files_pass    = Channel.fromPath("${params.run_dir}/pod5_pass/*.pod5")
    ch_pod5_files_fail    = Channel.fromPath("${params.run_dir}/pod5_fail/*.pod5")
    ch_pod5_files_skipped = Channel.fromPath("${params.run_dir}/pod5_skipped/*.pod5")
    ch_pod5_files         = ch_pod5_files_pass.mix(ch_pod5_files_fail).mix(ch_pod5_files_skipped).mix(ch_pod5_files)
    ch_collected_pod5     = ch_pod5_files.collect().ifEmpty([])

    //
    // CHANNEL: Put all pod5 generated files and their corresponding sample IDs into a single channel 
    //
    ch_pod5_files = ch_pod5_files
        .collate(params.dorado_batch_num)
        .map{ [[ id: it[0].simpleName.substring(0, 26) ], it ] }

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
    ch_meta = ch_meta.map{
        it.run_id = runid 
        it.id = it.sample_id
        it.remove("sample_id")
        it
    }

    //
    // CHANNEL: Collect barcode names
    //
    ch_barcodes = ch_meta.map{ it.barcode }.toSortedList()

    if (params.run_basecaller) {
        //
        // MODULE: Generate a bam file using pod5 files and any supplied bam to resume from
        //
        DORADO_BASECALLER (
            ch_pod5_files,
            params.bam ? ch_bam.map{it[1]} : [],
            dorado_model,
            params.dorado_bc_kit ?: []
        )
        ch_versions = ch_versions.mix(DORADO_BASECALLER.out.versions)
        ch_bam      = DORADO_BASECALLER.out.bam

        //
        // CHANNEL: Create basecalling merge channels
        //
        ch_bc_merge = ch_bam
            .collect{ it[1] }
            .branch {
                tomerge: it.size() > 1
                    return [[ id: it[0].simpleName.substring(0, 26) ], it ]
                pass: true
                    return [[ id: it[0].simpleName.substring(0, 26) ], it ]
            }

        //
        // MODULE: Merged basecalled bams if required
        //
        MERGE_BASECALLING (
            ch_bc_merge.tomerge,
            [[],[]],
            [[],[]]
        )
        ch_versions = ch_versions.mix(MERGE_BASECALLING.out.versions)
        ch_bam      = MERGE_BASECALLING.out.bam.mix(ch_bc_merge.pass)
    }

    if (params.run_demux) {
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
        // CHANNEL: Merge metadata to the demultiplexed fastq file
        //
        ch_demux_fastq = ch_meta
            .map { [it.barcode, it] }
            .join( ch_demux_fastq.map{it[1]}.flatten().map{ [ it.simpleName, it ] } )
            .map { [ it[1], it[2] ] }

        //
        // CHANNEL: Merge metadata to the demultiplexed bam file
        //
        ch_demux_bam = ch_meta
            .map { [it.barcode, it] }
            .join( ch_demux_bam.map{it[1]}.flatten().map{ [ it.simpleName, it ] } )
            .map { [ it[1], it[2] ] }
    }

    if (params.run_filtering) {
        if(params.emit_bam) {
            //
            // MODULE: Filter on read quality for bam files
            //
            FILTER_READ_Q (
                ch_demux_bam,
                [[], []],
                []
            )
            ch_versions  = ch_versions.mix(FILTER_READ_Q.out.versions)
            ch_demux_bam = FILTER_READ_Q.out.bam
        }
        else {
            //
            // MODULE: Filter on read quality for fastq files
            //
            CHOPPER (
                ch_demux_fastq
            )
            ch_versions             = ch_versions.mix(CHOPPER.out.versions)
            ch_demux_filtered_fastq = CHOPPER.out.fastq
        }
    }

    if(params.run_qc) {
        ch_fastqc_zip     = Channel.empty()
        ch_grouped_fastqc = Channel.empty()
        if(!params.emit_bam) {
            //
            // MODULE: Run fastqc on all fastq files if they exist
            //
            FASTQC (
                ch_demux_fastq
            )
            ch_versions   = ch_versions.mix(FASTQC.out.versions)
            ch_fastqc_zip = FASTQC.out.zip

            //
            // CHANNEL: Group fastqc by run_id,group,user,project
            //
            ch_grouped_fastqc = ch_fastqc_zip
                .map{ [ it[0].run_id, it[0].group, it[0].user, it[0].project_id, it ] }
                .groupTuple(by: [0, 1, 2, 3])
                .map {
                    def files = it[4].flatten().findAll { item -> !(item instanceof Map) }
                    [ [ id:it[1]+"_"+it[2]+"_"+it[3], run_id:it[0], group: it[1], user:it[2], project_id:it[3]], files ]
                }
        }

        //
        // MODULE: Run toulligqc on all samples
        //
        TOULLIGQC_ALL (
            ch_bam,
            ch_collected_pod5,
            ch_barcodes
        )
        ch_versions = ch_versions.mix(TOULLIGQC_ALL.out.versions)

        //
        // CHANNEL: Select seq data
        //
        if(params.emit_bam) {
            ch_seq = ch_demux_bam
        } else {
            ch_seq = ch_demux_fastq
        }

        //
        // CHANNEL: Extract metadata and bind to single bam file, split out barcodes from unclassified
        //
        ch_meta_bam = ch_seq
            .map { it[0] }
            .combine(ch_bam)
            .map { [ it[0], it[2] ] }
            .branch {
                bc: it[0].barcode != 'unclassified'
                uc: it[0].barcode == 'unclassified'
            }

        FILTER_BC (
            ch_meta_bam.bc,
            [[], []],
            []
        )
        ch_versions = ch_versions.mix(FILTER_BC.out.versions)
        ch_filt_bam = FILTER_BC.out.bam
        ch_filt_bam = ch_filt_bam.mix(ch_meta_bam.uc)

        //
        // CHANNEL: Group bam by run_id,group,user,project
        //
        ch_grouped_bam = ch_filt_bam
            .map{ [ it[0].run_id, it[0].group, it[0].user, it[0].project_id, it] }
            .groupTuple(by: [0, 1, 2, 3])
            .map {
                def files = it[4].flatten().findAll { item -> !(item instanceof Map) }
                [ [ id:it[1]+"_"+it[2]+"_"+it[3], run_id:it[0], group: it[1], user:it[2], project_id:it[3]], files ]
            }

        //
        // MODULE: Merge bams together
        //
        MERGE_GROUPS (
            ch_grouped_bam,
            [[], []],
            [[], []]
        )
        ch_versions = ch_versions.mix(MERGE_GROUPS.out.versions)
        ch_merged_bam = MERGE_GROUPS.out.bam

        //
        // MODULE: Create summary of merged bam
        //
        DORADO_SUMMARY (
            ch_merged_bam
        )
        ch_versions = ch_versions.mix(DORADO_SUMMARY.out.versions)

        //
        // CHANNEL: Filter out all bams with nothing in them
        //
        ch_merged_bam = ch_merged_bam
            .filter { row ->
                file(row[1]).size() >= 800
            }

        //
        // MODULE: Run toulligqc on merged bam
        //
        TOULLIGQC_GROUPED (
            ch_merged_bam,
            ch_collected_pod5,
            ch_barcodes
        )
        ch_versions = ch_versions.mix(TOULLIGQC_ALL.out.versions)

        //
        // CHANNEL: Group meta data by run_id,group,user,project and write to file
        //
        ch_grouped_meta = ch_filt_bam
            .map{
                def newMap = it[0].collectEntries { k, v -> [(k): v] }
                newMap.sample_id = newMap.id
                newMap.remove("id")
                newMap = [("sample_id"): newMap.remove("sample_id")] + newMap
                [newMap, it[1]]
            }
            .map{ [ it[0].run_id, it[0].group, it[0].user, it[0].project_id, it[0]] }
            .groupTuple(by: [0, 1, 2, 3])
            .map {
                def meta = it[4]
                def path = "${params.outdir}/grouped/${it[1]}/${it[2]}/asf/${it[3]}/${it[0]}/samplesheet.csv"
                dump_meta(meta, path)
            }

        //
        // MODULE: Collect software versions
        //
        CUSTOM_DUMPSOFTWAREVERSIONS (
            ch_versions.unique().collectFile()
        )

        // 
        // MODULE: MULTIQC
        // 
        workflow_summary = multiqc_summary(workflow, params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_unique_yml.collect())

        ch_multiqc_files_all = Channel.empty()
        ch_multiqc_files_all = ch_multiqc_files_all.mix(ch_multiqc_files)
        ch_multiqc_files_all = ch_multiqc_files_all.mix(ch_fastqc_zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files_all = ch_multiqc_files_all.collect().map{ [[id:"all"], it] }

        ch_multiqc_files_grouped = ch_grouped_fastqc
            .map { [ it[0].id, it[0], it ] }
            .map { [ it[1], [ it[2][1].flatten() ].flatten() ] }
            .combine(ch_multiqc_files.collect())
            .map { [ it[0], [it[1], it[2]].flatten() ] }

        MULTIQC_ALL (
            ch_multiqc_files_all,
            ch_multiqc_config,
            [],
            []
        )
        multiqc_report_all = MULTIQC_ALL.out.report.toList()

        MULTIQC_GROUPED (
            ch_multiqc_files_grouped,
            ch_multiqc_config,
            [],
            []
        )
        multiqc_report_grouped = MULTIQC_GROUPED.out.report.toList()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    dump_parameters(workflow, params)
    if(workflow.success) {
        workflow_complete_summary(workflow, "${params.outdir}/pipeline_info/workflow_complete.txt")
    }

    if (params.hook_url) {
        im_notification(workflow, params, projectDir, runid, summary_params, log)
    }

    // if (params.email || params.email_on_fail) {
    //     NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report, pass_mapped_reads, pass_trimmed_reads, pass_strand_check)
    // }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
