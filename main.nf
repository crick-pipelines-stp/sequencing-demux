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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DORADO_BASECALLER                } from './modules/local/dorado/basecaller/main'
include { DORADO_DEMUX                     } from './modules/local/dorado/demux/main'
include { SAMTOOLS_VIEW                    } from './modules/nf-core/samtools/view/main'
include { CHOPPER                          } from './modules/nf-core/chopper/main'
include { SEQKIT_SEQ                       } from './modules/nf-core/seqkit/seq/main' 
include { LINUX_COMMAND as FILTER_QC_FASTQ } from './modules/local/linux/command'
include { CAT_CAT as CAT_READ_IDS          } from './modules/nf-core/cat/cat/main'  
include { FASTQC                           } from './modules/nf-core/fastqc/main'
include { NANOPLOT as NANOPLOT_ALL         } from './modules/nf-core/nanoplot/main'
include { NANOPLOT as NANOPLOT_GROUPED     } from './modules/nf-core/nanoplot/main'
include { PYCOQC as PYCOQC_ALL             } from './modules/nf-core/pycoqc/main'
include { PYCOQC as PYCOQC_GROUPED         } from './modules/nf-core/pycoqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS      } from './modules/local/custom_dumpsoftwareversions'
include { MULTIQC as MULTIQC_ALL           } from './modules/local/multiqc/main'
include { MULTIQC as MULTIQC_GROUPED       } from './modules/local/multiqc/main'

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

    if (params.run_basecaller && !params.bam) {
        //
        // MODULE: Generate a bam file using the Dorado basecaller unless a bam file was already present as an input
        //
        DORADO_BASECALLER (
            ch_pod5_files,
            dorado_model
        )
        ch_versions = ch_versions.mix(DORADO_BASECALLER.out.versions)
        ch_bam      = DORADO_BASECALLER.out.bam
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
            SAMTOOLS_VIEW (
                ch_demux_bam,
                [],
                []
            )
            ch_versions  = ch_versions.mix(SAMTOOLS_VIEW.out.versions)
            ch_demux_bam = SAMTOOLS_VIEW.out.bam
        }
        else {
            //
            // MODULE: Filter on read quality for fastq files
            //
            CHOPPER (
                ch_demux_fastq
            )
            ch_versions    = ch_versions.mix(CHOPPER.out.versions)
            ch_demux_fastq = CHOPPER.out.fastq

        }
    }

    if(params.run_qc) {

        //
        // MODULE: Extract read ids from fastq file
        //
        SEQKIT_SEQ (
            ch_demux_fastq
        )
        ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions)
        ch_read_ids = SEQKIT_SEQ.out.fastx

        //
        // CHANNEL: Group read ids by run_id,group,user,project
        //
        ch_grouped_read_ids = ch_read_ids
            .map{ [ it[0].run_id, it[0].group, it[0].user, it[0].project_id, it ] }
            .groupTuple(by: [0, 1, 2, 3])
            .map {
                def files = it[4].flatten().findAll { item -> !(item instanceof Map) }
                [ [ id:it[1]+"_"+it[2]+"_"+it[3], run_id:it[0], group: it[1], user:it[2], project_id:it[3]], files ]
            }

        //
        // MODULE: Merge read ids from the same group
        //
        CAT_READ_IDS (
            ch_grouped_read_ids
        )
        ch_versions         = ch_versions.mix(CAT_READ_IDS.out.versions)
        ch_grouped_read_ids = CAT_READ_IDS.out.file_out

        //
        // MODULE: Filter sequencing summary file based on reads from groups
        //
        FILTER_QC_FASTQ (
            ch_grouped_read_ids,
            ch_sequencing_summary.collect()
        )
        ch_sequencing_summary_grouped = FILTER_QC_FASTQ.out.file.map{ [ it[0], it[1][1] ] }

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
        // MODULE: Run Nanoplot on all samples
        //
        NANOPLOT_ALL (
            ch_sequencing_summary
        )
        ch_versions = ch_versions.mix(NANOPLOT_ALL.out.versions)

        //
        // MODULE: Run Nanoplot on grouped samples
        //
        NANOPLOT_GROUPED (
            ch_sequencing_summary_grouped
        )
        ch_versions = ch_versions.mix(NANOPLOT_GROUPED.out.versions)

        //
        // MODULE: Run pycoqc on all samples
        //
        PYCOQC_ALL (
            ch_sequencing_summary
        )
        ch_versions = ch_versions.mix(PYCOQC_ALL.out.versions)

        //
        // MODULE: Run pycoqc on grouped samples
        //
        PYCOQC_GROUPED (
            ch_sequencing_summary_grouped
        )
        ch_versions = ch_versions.mix(PYCOQC_GROUPED.out.versions)

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
        ch_multiqc_files_all = ch_multiqc_files_all.mix(PYCOQC_ALL.out.json.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files_all = ch_multiqc_files_all.mix(NANOPLOT_ALL.out.txt.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files_all = ch_multiqc_files_all.collect().map{ [[id:"all"], it] }

        ch_multiqc_files_grouped = ch_grouped_fastqc
            .map { [ it[0].id, it[0], it ] }
            .join ( PYCOQC_GROUPED.out.json.map{ [ it[0].id, it[1] ] } )
            .join ( NANOPLOT_GROUPED.out.txt.map{ [ it[0].id, it[1] ] } )
            .map { [ it[1], [ it[2][1].flatten(), it[3], it[4] ].flatten() ] }
            .combine(ch_multiqc_files.collect())
            .map { [ it[0], [it[1], it[2], it[3], it[4]].flatten() ] }

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
