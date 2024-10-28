process CIRCOS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://docker.io/rlinder02/deciphergvizmsadatatable:v0.0.1':
        'docker.io/rlinder02/deciphergvizmsadatatable:v0.0.1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    generate_circos_plots.R $fasta $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //g')
    END_VERSIONS
    """
}
