/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { QCALIGN                } from '../subworkflows/local/qcalign'
include { FILTERCLUSTERCONSENSUS } from '../subworkflows/local/filterclusterconsensus'

include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_hifiampliconanalysis_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HIFIAMPLICONANALYSIS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_ref = ch_samplesheet.map { meta, fastq, fasta, primer1, primer2, bed, introns -> [meta, fasta] }
    ch_introns = ch_samplesheet.map { meta, fastq, fasta, primer1, primer2, bed, introns -> [meta, introns] }
    ch_fastq = ch_samplesheet.map { meta, fastq, fasta, primer1, primer2, bed, introns -> [meta, fastq] }
    // only process fastq.gz or fq.gz files (reads); fasta files get processed later
    ch_fastq_only = ch_fastq.map { meta, file -> 
                    def fileType = file.name.toString().split('/').last().split('\\.').last()
                    if (fileType != "fasta") {
                        return tuple(meta, file)
                    }
                 }

    //
    // SUBWORKFLOW: Align HiFi reads to gene-specific genome and run QC on raw and aligned reads
    //
    QCALIGN ( ch_fastq_only,
                ch_ref,
                ch_introns
            )

    ch_multiqc_files = ch_multiqc_files.mix(QCALIGN.out.nanostats_report.map {it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(QCALIGN.out.bam_qc.map {it[1]})
    ch_versions = ch_versions.mix(QCALIGN.out.versions)
    aligned_fasta = QCALIGN.out.aligned_fasta
    
    //
    // SUBWORKFLOW: Filter for aligned reads with both primers sequences present, then cluster based on sequence similarity, then identify a single consensus sequence for each cluster 
    //
    FILTERCLUSTERCONSENSUS ( ch_samplesheet,
                             aligned_fasta
                           )

    ch_versions = ch_versions.mix(FILTERCLUSTERCONSENSUS.out.versions)
    
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
