include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check/main'

workflow SAMPLESHEET_PARSE {
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
        .splitCsv (header:true, sep:",")
        .map {
            it.group = it.group.replaceAll(" ", "_").toLowerCase()
            it.user = it.user.replaceAll(" ", "_").toLowerCase()
            if(params.dorado_bc_parse_pos != null && params.dorado_bc_kit != null) {
                it.barcode = "barcode" + it.barcode.substring(params.dorado_bc_parse_pos, (params.dorado_bc_parse_pos + 2)
            }
            if(params.dorado_bc_kit != null) {
                it.barcode = params.dorado_bc_kit + "_" + it.barcode
            }
            it
        }

    emit:
    meta // channel: [ val(meta) ]
    versions = SAMPLESHEET_CHECK.out.versions
}
