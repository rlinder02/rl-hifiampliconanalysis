include { CONVERTTOFASTA } from '../../modules/local/converttofasta'
include { FCS_FCSADAPTOR } from '../../modules/nf-core/fcs/fcsadaptor/main'
include { NANOPLOT       } from '../../modules/nf-core/nanoplot/main'

workflow QCALIGN {

    take:
    ch_samplesheet 

    main:

    ch_fastq = ch_samplesheet.map { meta, fastq, fasta, primer1, primer2, bed -> [meta, fastq] }
    ch_versions = Channel.empty()

    CONVERTTOFASTA ( ch_fastq )
    ch_versions = ch_versions.mix(CONVERTTOFASTA.out.versions.first())

    FCS_FCSADAPTOR ( CONVERTTOFASTA.out.fasta )
    ch_versions = ch_versions.mix(FCS_FCSADAPTOR.out.versions.first())

    NANOPLOT ( ch_fastq )
    ch_versions = ch_versions.mix(NANOPLOT.out.versions.first())


    emit:
    vec_report            = FCS_FCSADAPTOR.out.adaptor_report      // channel: [ val(meta), [ report ] ]
    nanostats_report      = NANOPLOT.out.txt                       // channel: [ val(meta), [ txt ] ]
    nanostats_results     = NANOPLOT.out.nanoplot_results          // channel: [ val(meta), [ results_dir ] ]
    versions              = ch_versions                            // channel: [ versions.yml ]
}

