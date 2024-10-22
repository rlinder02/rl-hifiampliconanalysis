process CALLCONSENSUS {
    tag "$meta.id"
    label 'process_medium'
    debug true

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
    file_name=\$(basename $bam .bam)
    echo \$file_name
    cluster_id=\${file_name##*_}
    echo \$cluster_id

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
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
