process DORADO_DEMUX {
    tag "$meta.id"
    label 'process_high'

    container 'thecrick/pipetech_dorado:0.6.0-linux-x64'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam") , emit: bam
    path  "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    dorado demux \\
        --output-dir ./ \\
        --no-classify \\
        $args \\
        $bam 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(dorado --version 2>&1)
    END_VERSIONS
    """
}
