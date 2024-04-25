include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check/main'

workflow CHECK_SAMPLESHEET {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    //
    // MODULE: Run python samplesheet checker
    //
    SAMPLESHEET_CHECK ( 
        samplesheet 
    )

    //
    // CHANNEL: Load csv and split columns into metadata
    //
    meta = SAMPLESHEET_CHECK.out.csv
        .splitCsv ( header:true, sep:"," )

    emit:
    meta // channel: [ val(meta) ]
    versions = SAMPLESHEET_CHECK.out.versions
}