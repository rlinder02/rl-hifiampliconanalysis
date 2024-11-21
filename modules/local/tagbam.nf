process TAGBAM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_biopython:1938b0dd1fb1aab7':
        'community.wave.seqera.io/library/pysam_biopython:1938b0dd1fb1aab7' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*modified.bam"), emit: modified_bam
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tag_clusters.py $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
