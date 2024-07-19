process DORADO_BASECALLER {
    tag "$meta.id"
    label 'process_high'

    container 'docker.io/thecrick/pipetech_dorado:0.7.1-linux-x64'

    input:
    tuple val(meta), path("pod5s/*")
    val(model)
    val(bc_kit)

    output:
    tuple val(meta), path("*.bam") , emit: bam
    path  "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def bc_kit_arg = bc_kit ? "--kit-name ${bc_kit}"  : ''

    """
    dorado basecaller \\
        $model \\
        pod5s/ \\
        $bc_kit_arg \\
        $args \\
        > ${prefix}.bam


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(dorado --version 2>&1)
    END_VERSIONS
    """
}
