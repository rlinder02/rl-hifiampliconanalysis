process SPECIESPLOTS {
    tag "$meta"
    label 'process_low'
    debug true

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://docker.io/rlinder02/ggtranscript:v0.0.5':
        'docker.io/rlinder02/ggtranscript:v0.0.5' }"

    input:
    tuple val(meta), path(vcfs), path(bed), path(bounds), path(total_reads), path(orfs)

    output:
    tuple val(meta), path("*.png"), emit: png
    tuple val(meta), path("*.svg"), emit: svg
    tuple val(meta), path("*.bed"), emit: struct_bed, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    """
    vcf_list=\$(echo $vcfs | sed 's/ /\\n/g' | grep 'Preprocessed' | sort > vcf_fofn.txt)
    vcf_list2=\$(echo $vcfs | sed 's/ /\\n/g' | grep -v 'Preprocessed' | sort >> vcf_fofn.txt)
    orf_list=\$(echo $orfs | sed 's/ /\\n/g' | grep 'Preprocessed' | sort > orf_fofn.txt)
    orf_list2=\$(echo $orfs | sed 's/ /\\n/g' | grep -v 'Preprocessed' | sort >> orf_fofn.txt)
    total_reads_list=\$(echo $total_reads | sed 's/ /\\n/g' | grep 'Preprocessed' | sort > total_reads_fofn.txt)
    total_reads_list2=\$(echo $total_reads | sed 's/ /\\n/g' | grep -v 'Preprocessed' | sort >> total_reads_fofn.txt)

    generate_linear_plots.R vcf_fofn.txt $bed $bounds total_reads_fofn.txt $meta orf_fofn.txt

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
