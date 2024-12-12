process CALLCONSENSUSPP {
    tag "$meta.id"
    label 'process_medium'
    debug true

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bcftools_crossmap_minimap2_orfipy_pruned:a10b34d7e4795d31':
        'community.wave.seqera.io/library/bcftools_crossmap_minimap2_orfipy_pruned:a10b34d7e4795d31' }"
    containerOptions = "--user root"

    input:
    tuple val(meta), path(bam), path(ref)

    output:
    path("*modified.vcf.gz")                                     , emit: vcf        , optional: true
    path("*.fasta")                                              , emit: con_fasta  , optional: true
    path("orfipy/*.bed")                                         , emit: orf_bed    , optional: true
    path("*.bed")                                                , emit: lifted_bed , optional: true
    tuple val(meta), path("*.txt")                               , emit: txt        , optional: true
    path("*chain")                                               , emit: chain      , optional: true
    path("*.paf")                                                , emit: paf        , optional: true
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    file_name=\$(basename $bam .bam)
    cluster_id=\${file_name##*_}
    echo \$cluster_id > cluster_id
    samtools \\
        view \\
        $bam \\
    | \\
    cut -f1 | sort | uniq | wc -l > ${prefix}_\${cluster_id}_aligned_reads.txt

    bcftools \\
        mpileup \\
        $args \\
        --threads $task.cpus \\
        --max-depth 1000000 \\
        -L 1000000 \\
        -X pacbio-ccs \\
        -Ou \\
        -f $ref \\
        $bam \\
    | \\
    bcftools \\
        call \\
        -m \\
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
        -i 'QUAL >= 20' \\
        -Oz \\
        -s FAIL \\
        --threads $task.cpus \\
        -o ${prefix}_\${cluster_id}_modified.vcf.gz

    tabix -p vcf ${prefix}_\${cluster_id}_modified.vcf.gz
    
    if [[ \$(zcat ${prefix}_\${cluster_id}_modified.vcf.gz | grep 'PASS' | grep -c 'AC=') -gt 0 || \$(cat $ref | grep -v '>' | tr -d '\\n' | wc -m) -ne \$(zcat ${prefix}_\${cluster_id}_modified.vcf.gz | grep -v '#' | grep -vc 'DP=0') ]]
    then
        bcftools \\
            consensus \\
            -a '*' \\
            -i 'QUAL >= 20' \\
            -o ${prefix}_\${cluster_id}.fasta \\
            -f $ref \\
            ${prefix}_\${cluster_id}_modified.vcf.gz
        
        cat ${prefix}_\${cluster_id}.fasta | tr -d '*' | tr -d '-' | tr -d '\\n' | sed 's/_cDNA/_cDNA\\n/g' > ${prefix}_\${cluster_id}_modified.fasta

        orfipy \\
            ${prefix}_\${cluster_id}_modified.fasta \\
            --bed ${prefix}_\${cluster_id}.bed \\
            --outdir orfipy \\
            --procs $task.cpus

        minimap2 \\
            -cx asm20 \\
            --cs \\
            ${prefix}_\${cluster_id}_modified.fasta \\
            ${prefix}_\${cluster_id}.fasta \\
            -o ${prefix}_\${cluster_id}_modified.paf
        
    else
        orfipy \\
            $ref \\
            --bed ${prefix}_\${cluster_id}.bed \\
            --outdir orfipy \\
            --procs $task.cpus
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
