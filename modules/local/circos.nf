process CIRCOS {
    tag "$meta"
    label 'process_low'
    debug true

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-complexheatmap_r-circlize_r-data.table_r-gridbase_pruned:8987b0b2d3ee3ff2':
        'community.wave.seqera.io/library/bioconductor-complexheatmap_r-circlize_r-data.table_r-gridbase_pruned:8987b0b2d3ee3ff2' }"

    input:
    tuple val(meta), path(vcfs), path(bed), path(bounds), path(total_reads)

    output:
    tuple val(meta), path("*.png"), emit: png
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def vcfList = vcfs.join(' ')
    """
    vcf_list=\$(echo $vcfs | sed 's/ /\n/g')
    echo \$vcf_list
    
    generate_circos_plots.R \$vcf_list $bed $bounds $total_reads $meta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //g')
    END_VERSIONS
    """
}
