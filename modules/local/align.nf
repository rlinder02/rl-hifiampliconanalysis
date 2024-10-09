process ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/minimap2_samtools:7e38c0cfb1291cfb':
        'community.wave.seqera.io/library/minimap2_samtools:7e38c0cfb1291cfb' }"

    input:
    tuple val(meta), path(fasta), path(ref)

    output:
    tuple val(meta), path("*.sorted.fasta"), emit: sorted_fasta
    tuple val(meta), path("*.sorted.bam")  , emit: sorted_bam
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    minimap2 \\
        --eqx \\
        --sam-hit-only \\
        -t $task.cpus \\
        -Yax splice:hq \\ 
        $args \\
        $ref \\
        $fasta \\
    | \\
    samtools \\
        sort \\
        $args \\
        -m4G \\
        -@ $task.cpus \\
        -O BAM \\
        -o ${prefix}.sorted.bam
    
    samtools \\
        fasta \\
        -@ $task.cpus \\
        ${prefix}.sorted.bam > ${prefix}.sorted.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        align: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        align: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
