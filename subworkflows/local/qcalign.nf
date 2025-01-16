include { CONVERTTOFASTA } from '../../modules/local/converttofasta'
include { FCS_FCSADAPTOR } from '../../modules/nf-core/fcs/fcsadaptor/main'
include { NANOPLOT       } from '../../modules/nf-core/nanoplot/main'
include { ALIGNINTRONS   } from '../../modules/local/alignintrons'
include { ALIGN          } from '../../modules/local/align'
include { QUALIMAP_BAMQC } from '../../modules/nf-core/qualimap/bamqc/main'

workflow QCALIGN {

    take:
    ch_fastq_only
    ch_ref
    ch_introns

    main:
    ch_versions = Channel.empty()

    CONVERTTOFASTA ( ch_fastq_only )
    ch_versions = ch_versions.mix(CONVERTTOFASTA.out.versions.first())

    // need to change FCS_FCSADAPTOR so if no contaminating sequences found, still outputs a fa.gz file for alignment step
    FCS_FCSADAPTOR ( CONVERTTOFASTA.out.fasta )
    ch_versions = ch_versions.mix(FCS_FCSADAPTOR.out.versions.first())

    NANOPLOT ( ch_fastq_only )
    ch_versions = ch_versions.mix(NANOPLOT.out.versions.first())

    ch_intron_fasta = FCS_FCSADAPTOR.out.cleaned_assembly.combine(ch_introns, by:0)
    
    ALIGNINTRONS ( ch_intron_fasta )
    ch_versions = ch_versions.mix(ALIGNINTRONS.out.versions.first())

    ch_intronless_ref = ALIGNINTRONS.out.sorted_fasta.combine(ch_ref, by:0)

    ALIGN ( ch_intronless_ref )
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