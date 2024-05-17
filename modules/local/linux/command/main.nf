process LINUX_COMMAND {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(input2)

    output:
    tuple val(meta), path("*.cmd.*"), emit: file
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def ext      = task.ext.ext ?: 'txt'
    def cmd      = task.ext.cmd ?: 'echo "NO-ARGS"'
    def pre_cmd  = task.ext.pre_cmd ?: ''
    def post_cmd = task.ext.post_cmd ?: ''

    """
    $pre_cmd
    cat $input | $cmd > ${prefix}.cmd.${ext}
    $post_cmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        linux: NOVERSION
    END_VERSIONS
    """
}
