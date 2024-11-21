process FINDPRIMERS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_biopython:1938b0dd1fb1aab7':
        'community.wave.seqera.io/library/pysam_biopython:1938b0dd1fb1aab7' }"

    input:
    tuple val(meta), path(fasta), val(primer1), val(primer2)

    output:
    tuple val(meta), path("*_filtered.fasta"), emit: filtered_fasta
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    find_primers.py $fasta -p1 $primer1 -p2 $primer2 -t $task.cpus -m 2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
