process CIRCOS {
    tag "$meta"
    label 'process_low'
    debug true

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-complexheatmap_r-circlize_r-data.table_r-gridbase_pruned:8987b0b2d3ee3ff2':
        'community.wave.seqera.io/library/bioconductor-complexheatmap_r-circlize_r-data.table_r-gridbase_pruned:8987b0b2d3ee3ff2' }"

    input:
    tuple val(meta), path(vcfs), path(bed), path(bounds), path(total_reads), path(orfs)

    output:
    tuple val(meta), path("*.png"), emit: png
    tuple val(meta), path("*.bed"), emit: struct_bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    """
    vcf_list=\$(echo $vcfs | sed 's/ /\\n/g' | sort > vcf_fofn.txt)
    orf_list=\$(echo $orfs | sed 's/ /\\n/g' | sort > orf_fofn.txt)

    generate_circos_plots.R vcf_fofn.txt $bed $bounds $total_reads $meta orf_fofn.txt

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
