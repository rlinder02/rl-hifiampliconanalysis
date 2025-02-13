process ALIGNINTRONS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/minimap2_samtools:7e38c0cfb1291cfb':
        'community.wave.seqera.io/library/minimap2_samtools:7e38c0cfb1291cfb' }"

    input:
    tuple val(meta), path(fasta), path(introns)

    output:
    tuple val(meta), path("*.intronless.sorted.fasta"), emit: sorted_fasta
    tuple val(meta), path("*.intronless.sorted.bam")  , emit: sorted_bam
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    minimap2 \\
        --eqx \\
        -t $task.cpus \\
        -Y \\
        -a \\
        -x splice:hq \\
        $introns \\
        $fasta \\
    | \\
    samtools \\
        view \\
        -f 4 \\
    | \\
    samtools \\
        sort \\
        -m4G \\
        -@ $task.cpus \\
        -O BAM \\
        -o ${prefix}.intronless.sorted.bam
    
    samtools \\
        fasta \\
        -@ $task.cpus \\
        ${prefix}.intronless.sorted.bam > ${prefix}.intronless.sorted.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
