process CONVERTTOFASTA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/seqtk:1.4--2691a67b9f188d94':
        'community.wave.seqera.io/library/seqtk:1.4--2691a67b9f188d94' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def SEQTK_VERSION = '1.4-r122' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    seqtk seq \\
        -a $fastq \\
        $args \\
        > ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: $SEQTK_VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def SEQTK_VERSION = '1.4-r122' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: $SEQTK_VERSION
    END_VERSIONS
    """
}
