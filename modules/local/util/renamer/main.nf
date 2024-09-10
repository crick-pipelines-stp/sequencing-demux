process UTIL_RENAMER {
    tag "$meta.id"
    executor 'local'

    input:
    tuple val(meta), path(file)
    val(rename)

    output:
    tuple val(meta), path("*.*", includeInputs: true), emit: output

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ext = task.ext.ext ?: ".noext"
    def command = rename ? "mv ${file} ${prefix}.${ext}" : ''

    """
    $command
    """
}
