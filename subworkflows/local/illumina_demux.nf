/*
Subworkflow for demultiplexing illumina data
*/

workflow ILLUMINA_DEMULTIPLEX {

    take:
    val_samplesheet   // The string path of the samplesheet to parse for metadata

    main:

    // Init
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // ch_file_1 = Channel.of(file("/Users/cheshic/dev/test_data/illumina/bcl_convert/Reports_1/Adapter_Metrics.csv", checkIfExists: true))
    // ch_file_2 = Channel.of(file("/Users/cheshic/dev/test_data/illumina/bcl_convert/Reports_1/Demultiplex_Stats.csv", checkIfExists: true))
    // ch_file_3 = Channel.of(file("/Users/cheshic/dev/test_data/illumina/bcl_convert/Reports_1/Quality_Metrics.csv", checkIfExists: true))
    // ch_file_4 = Channel.of(file("/Users/cheshic/dev/test_data/illumina/bcl_convert/Reports_1/RunInfo.xml", checkIfExists: true))
    // ch_file_5 = Channel.of(file("/Users/cheshic/dev/test_data/illumina/bcl_convert/Reports_1/Top_Unknown_Barcodes.csv", checkIfExists: true))

    // ch_file_6 = Channel.of(file("/Users/cheshic/dev/test_data/illumina/bcl_convert/Reports_2/Adapter_Metrics.csv", checkIfExists: true))
    // ch_file_7 = Channel.of(file("/Users/cheshic/dev/test_data/illumina/bcl_convert/Reports_2/Demultiplex_Stats.csv", checkIfExists: true))
    // ch_file_8 = Channel.of(file("/Users/cheshic/dev/test_data/illumina/bcl_convert/Reports_2/Quality_Metrics.csv", checkIfExists: true))
    // ch_file_9 = Channel.of(file("/Users/cheshic/dev/test_data/illumina/bcl_convert/Reports_2/RunInfo.xml", checkIfExists: true))
    // ch_file_10 = Channel.of(file("/Users/cheshic/dev/test_data/illumina/bcl_convert/Reports_2/Top_Unknown_Barcodes.csv", checkIfExists: true))

    ch_file_1 = Channel.of(file("/Users/cheshic/dev/test_data/illumina/bcl_convert/Reports_1", checkIfExists: true))
    ch_file_2 = Channel.of(file("/Users/cheshic/dev/test_data/illumina/bcl_convert/Reports_2", checkIfExists: true))

    // ch_multiqc_files = ch_file_1
    //     .mix(ch_file_2)
    //     .mix(ch_file_3)
    //     .mix(ch_file_4)
    //     .mix(ch_file_5)
    //     .mix(ch_file_6)
    //     .mix(ch_file_7)
    //     .mix(ch_file_8)
    //     .mix(ch_file_9)
    //     .mix(ch_file_10)

    ch_multiqc_files = ch_file_1
        .mix(ch_file_2)
        // .mix(ch_file_3)
        // .mix(ch_file_4)
        // .mix(ch_file_5)
        // .mix(ch_file_6)
        // .mix(ch_file_7)
        // .mix(ch_file_8)
        // .mix(ch_file_9)
        // .mix(ch_file_10)

    emit:
    versions      = ch_versions
    multiqc_files = ch_multiqc_files
}