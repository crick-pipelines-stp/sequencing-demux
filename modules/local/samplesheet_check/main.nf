process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "docker.io/python:3.11"

    input:
    path samplesheet

    output:
    path '*.csv'        , emit: csv
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    sample       = samplesheet
    process_name = task.process
    output       = task.ext.output ?: 'samplesheet.valid.csv'
    template 'samplesheet_check.py'
}
