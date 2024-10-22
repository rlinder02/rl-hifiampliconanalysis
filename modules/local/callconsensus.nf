process CALLCONSENSUS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bcftools_minimap2_samtools:c2489975f9638f9b':
        'community.wave.seqera.io/library/bcftools_minimap2_samtools:c2489975f9638f9b' }"

    input:
    tuple val(meta), path(bam), path(ref)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools \\
        mpileup \\
        $args \\
        --threads $task.cpus \\
        --max-depth 100000 \\
        -L 100000 \\
        -X pacbio-ccs \\
        -Ou \\
        -f $ref \\
        $bam \\
    | \\
    bcftools \\
        call \\
        -mv \\
        -Oz \\
        -o ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        callconsensus: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        callconsensus: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
