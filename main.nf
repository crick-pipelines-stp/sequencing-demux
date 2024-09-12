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

include { params_summary_map        } from './modules/francis-crick-institute/util/logging/main'
include { summary_log               } from './modules/francis-crick-institute/util/logging/main'
include { multiqc_summary           } from './modules/francis-crick-institute/util/logging/main'
include { im_notification           } from './modules/francis-crick-institute/util/logging/main'
include { dump_parameters           } from './modules/francis-crick-institute/util/io/main'
include { dump_meta                 } from './modules/francis-crick-institute/util/io/main'
include { gen_id                    } from './modules/francis-crick-institute/util/io/main'
include { workflow_complete_summary } from './modules/francis-crick-institute/util/io/main'

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
def check_param_list = [
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
    params.bam,
    params.dorado_resume_bam
]
for (param in check_param_list) { if (param) { file(param, checkIfExists: true) } }

// Check mode of operation
def operation_modes = ['ont', 'illumina']
if(!(params.mode in operation_modes)) {
    exit 1, "Incorrect mode of operation '${params.mode}'. Please choose from ${operation_modes}"
}

// Extract run_id
def runid = file(params.run_dir).name

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_VIEW as FILTER_READ_Q      } from './modules/nf-core/samtools/view/main'
include { UTIL_RENAMER as OUTPUT_RAW_BAM      } from './modules/local/util/renamer/main'
include { TOULLIGQC as TOULLIGQC_ALL          } from './modules/local/toulligqc/main'
include { SAMTOOLS_MERGE as MERGE_GROUPS      } from './modules/nf-core/samtools/merge/main'
include { DORADO_SUMMARY                      } from './modules/francis-crick-institute/dorado/summary/main'
include { TOULLIGQC as TOULLIGQC_GROUPED      } from './modules/local/toulligqc/main'
include { SAMTOOLS_FASTQ as RAW_BAM_TO_FASTQ  } from './modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_FASTQ as FILT_BAM_TO_FASTQ } from './modules/nf-core/samtools/fastq/main'
include { FASTQC                              } from './modules/nf-core/fastqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS         } from './modules/local/custom_dumpsoftwareversions'
include { MULTIQC as MULTIQC_ALL              } from './modules/local/multiqc/main'
include { MULTIQC as MULTIQC_GROUPED          } from './modules/local/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ONT_DEMULTIPLEX } from './subworkflows/francis-crick-institute/ont/demultiplex/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    //
    // INIT:
    // 
    ch_versions   = Channel.empty()
    ch_fastqc_zip = Channel.empty()

    if(params.mode == "ont")
    {
        //
        // MODULE: Demultiplex nanopore data
        //
        ONT_DEMULTIPLEX(
            params.dorado_model,
            params.dorado_bc_kit,
            true,
            params.dorado_bc_parse_pos,
            params.dorado_append_bc,
            params.dorado_batch_num,
            params.run_dir,
            params.bam,
            params.dorado_resume_bam,
            true,
            params.samplesheet
        )
        ch_versions       = ch_versions.mix(ONT_DEMULTIPLEX.out.versions)
        ch_pod5           = ONT_DEMULTIPLEX.out.pod5
        ch_basecalled_bam = ONT_DEMULTIPLEX.out.bam
        ch_demux_bam      = ONT_DEMULTIPLEX.out.demux_bam

        //
        // CHANNEL: Create specialised channels from demux subworkflow
        //
        ch_collected_pod5 = ch_pod5.collect{it[1]}.ifEmpty([])
        ch_barcodes = ch_demux_bam.map{ it[0].barcode }.toSortedList()
        ch_demux_bam = ch_demux_bam.map{
            it[0].run_id = runid
            it
        }
        ch_demux_bam.collect{it[0]}.map {
            dump_meta(it, "${params.outdir}/pipeline_info/sample_meta.csv")
        }

        //
        // MODULE: Output the raw demux bams if required
        //
        if(params.output_raw && params.output_bam)
        {
            OUTPUT_RAW_BAM (
                ch_demux_bam,
                false
            )
        }

        //
        // MODULE: Filter on read quality files
        //
        FILTER_READ_Q (
            ch_demux_bam,
            [[], []],
            []
        )
        ch_versions     = ch_versions.mix(FILTER_READ_Q.out.versions)
        ch_filtered_bam = FILTER_READ_Q.out.bam

        //
        // MODULE: Run toulligqc on all samples
        //
        TOULLIGQC_ALL (
            ch_basecalled_bam,
            ch_collected_pod5,
            ch_barcodes
        )
        ch_versions = ch_versions.mix(TOULLIGQC_ALL.out.versions)

        //
        // CHANNEL: Group bam by run_id,group,user,project
        //
        ch_grouped_bam = ch_demux_bam
            .map{ [ it[0].run_id, it[0].group, it[0].user, it[0].project_id, it] }
            .groupTuple(by: [0, 1, 2, 3])
            .map {
                def files = it[4].flatten().findAll { item -> !(item instanceof Map) }
                [ [ id:it[1]+"_"+it[2]+"_"+it[3], run_id:it[0], group: it[1], user:it[2], project_id:it[3]], files ]
            }

        //
        // MODULE: Merge bams together by their grouping
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
        // CHANNEL: Group meta data by run_id,group,user,project and write samplesheet to file
        //
        ch_grouped_meta = ch_demux_bam
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
        // MODULE: Convert raw bams to fastq
        //
        RAW_BAM_TO_FASTQ (
            ch_demux_bam,
            false
        )
        ch_versions = ch_versions.mix(RAW_BAM_TO_FASTQ.out.versions)
        ch_raw_fastq = RAW_BAM_TO_FASTQ.out.reads

        //
        // MODULE: Convert filtered bams to fastq
        //
        FILT_BAM_TO_FASTQ (
            ch_filtered_bam,
            false
        )
        ch_versions = ch_versions.mix(FILT_BAM_TO_FASTQ.out.versions)
        ch_filt_fastq = FILT_BAM_TO_FASTQ.out.reads

        //
        // MODULE: Run fastqc on fastq files
        //
        FASTQC (
            ch_raw_fastq
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
