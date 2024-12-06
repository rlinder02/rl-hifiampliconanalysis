include { CONVERTTOFASTA } from '../../modules/local/converttofasta'
include { FCS_FCSADAPTOR } from '../../modules/nf-core/fcs/fcsadaptor/main'
include { NANOPLOT       } from '../../modules/nf-core/nanoplot/main'
include { ALIGN          } from '../../modules/local/align'
include { QUALIMAP_BAMQC } from '../../modules/nf-core/qualimap/bamqc/main'

workflow QCALIGN {

    take:
    ch_samplesheet 

    main:

    ch_fastq = ch_samplesheet.map { meta, fastq, fasta, primer1, primer2, bed -> [meta, fastq] }
    ch_ref = ch_samplesheet.map { meta, fastq, fasta, primer1, primer2, bed -> [meta, fasta] }
    ch_versions = Channel.empty()
    is_fastq = ch_fastq.map { meta, file -> file.toString().contains('q.gz') ? 1: 0 }
    // only process fastq.gz or fq.gz files (reads); fasta files get processed later
    if (is_fastq == 1) {
        CONVERTTOFASTA ( ch_fastq )
        ch_versions = ch_versions.mix(CONVERTTOFASTA.out.versions.first())

        // need to change FCS_FCSADAPTOR so if no contaminating sequences found, still outputs a fa.gz file for alignment step
        FCS_FCSADAPTOR ( CONVERTTOFASTA.out.fasta )
        ch_versions = ch_versions.mix(FCS_FCSADAPTOR.out.versions.first())

        NANOPLOT ( ch_fastq )
        ch_versions = ch_versions.mix(NANOPLOT.out.versions.first())

        ch_fastas = FCS_FCSADAPTOR.out.cleaned_assembly.combine(ch_ref, by:0)

        ALIGN ( ch_fastas )
        ch_versions = ch_versions.mix(ALIGN.out.versions.first())

        QUALIMAP_BAMQC ( ALIGN.out.sorted_bam )
        ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions.first())

        emit:
        vec_report            = FCS_FCSADAPTOR.out.adaptor_report      // channel: [ val(meta), [ report ] ]
        nanostats_report      = NANOPLOT.out.txt                       // channel: [ val(meta), [ txt ] ]
        nanostats_results     = NANOPLOT.out.nanoplot_results          // channel: [ val(meta), [ results_dir ] ]
        aligned_fasta         = ALIGN.out.sorted_fasta                 // channel: [ val(meta), [ fasta ] ]
        bam_qc                = QUALIMAP_BAMQC.out.results             // channel: [ val(meta), [ outdir ] ]
        versions              = ch_versions                            // channel: [ versions.yml ]
    }
}

