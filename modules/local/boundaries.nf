process BOUNDARIES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_biopython:1938b0dd1fb1aab7':
        'community.wave.seqera.io/library/pysam_biopython:1938b0dd1fb1aab7' }"

    input:
    tuple val(meta), path(ref), val(primer1), val(primer2)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    primer_boundaries.py $ref -p1 $primer1 -p2 $primer2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
