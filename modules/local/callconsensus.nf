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
    tuple val(meta.id), path("*modified.vcf.gz"), emit: vcf       , optional: true
    tuple val(meta.id), path("*.fasta")         , emit: con_fasta , optional: true
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    file_name=\$(basename $bam .bam)
    cluster_id=\${file_name##*_}

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
        -Ou \\
        --threads $task.cpus \\
    | \\
    bcftools \\
        norm \\
        -f $ref \\
        --threads $task.cpus \\
        -Ou \\
    | \\
    bcftools \\
        filter \\
        --IndelGap 5 \\
        -i 'QUAL >= 20 & INFO/DP >= 5' \\
        -Oz \\
        --threads $task.cpus \\
        -o ${prefix}_\${cluster_id}.vcf.gz

    bcftools +setGT \\
        ${prefix}_\${cluster_id}.vcf.gz -- \\
        -t q \\
        -i 'AD[:1]/DP>=0.8' \\
        -n 'c:1/1' \\
    | \\
    bcftools +setGT \\
        -o ${prefix}_\${cluster_id}_modified.vcf.gz \\
        -- \\
        -t q \\
        -i 'AD[:1]/DP<0.8' \\
        -n 'c:0/1'
    tabix -p vcf ${prefix}_\${cluster_id}_modified.vcf.gz
    
    if [ \$(zcat ${prefix}_\${cluster_id}_modified.vcf.gz | grep -vc '#') -gt 0 ]
    then
        bcftools \\
            consensus \\
            -o ${prefix}_\${cluster_id}.fasta \\
            -f $ref \\
            -H I \\
            ${prefix}_\${cluster_id}_modified.vcf.gz
    else 
        rm -f ${prefix}_\${cluster_id}_modified.vcf.gz
    fi

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
