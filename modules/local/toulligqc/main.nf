process TOULLIGQC {
    label 'process_low'
    tag "$meta.id"

    container "docker.io/genomicpariscentre/toulligqc:2.7"

    input:
    tuple val(meta), path(ontfile)

    output:
    tuple val(meta), path("*/*.data")              , emit: report_data
    tuple val(meta), path("*/*.html")              , emit: report_html, optional: true
    tuple val(meta), path("*/images/*.html")       , emit: plots_html
    tuple val(meta), path("*/images/plotly.min.js"), emit: plotly_js
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def input_file = ("$ontfile".endsWith(".fastq") || "$ontfile".endsWith(".fastq.gz") || "$ontfile".endsWith(".fq") || "$ontfile".endsWith(".fq.gz")) ? "--fastq ${ontfile}" :
        ("$ontfile".endsWith(".txt") || "$ontfile".endsWith(".txt.gz")) ? "--sequencing-summary-source ${ontfile}" :
        ("$ontfile".endsWith(".bam")) ? "--bam ${ontfile}" : ''

    """
    toulligqc \\
        $input_file \\
        --output-directory ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        toulligqc: \$(toulligqc --version)
    END_VERSIONS
    """
}
