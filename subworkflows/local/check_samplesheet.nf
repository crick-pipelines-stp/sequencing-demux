include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check/main'

workflow CHECK_SAMPLESHEET {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )

    meta = SAMPLESHEET_CHECK.out.csv
        .splitCsv ( header:true, sep:"," )

    emit:
    meta // channel: [ val(meta) ]
    versions = SAMPLESHEET_CHECK.out.versions
}
