include { FINDPRIMERS      } from '../../modules/local/findprimers'
include { CLUSTER          } from '../../modules/local/cluster'
include { ALIGNCLUSTERS    } from '../../modules/local/alignclusters'
include { SPLITBAM         } from '../../modules/local/splitbam'

workflow FILTERCLUSTERCONSENSUS {

    take:
    ch_samplesheet 
    ch_aligned_fasta   // channel: [ val(meta), [ fasta ] ]

    main:

    ch_primer1 = ch_samplesheet.map { meta, fastq, fasta, primer1, primer2, bed -> [meta, primer1] }
    ch_primer2 = ch_samplesheet.map { meta, fastq, fasta, primer1, primer2, bed -> [meta, primer2] }
    ch_ref = ch_samplesheet.map { meta, fastq, fasta, primer1, primer2, bed -> [meta, fasta] }
    ch_versions = Channel.empty()

    ch_fasta_primer1 = ch_aligned_fasta.combine(ch_primer1, by:0)
    ch_fasta_primer1_primer2 = ch_fasta_primer1.combine(ch_primer2, by:0)
    
    FINDPRIMERS ( ch_fasta_primer1_primer2 )
    ch_versions = ch_versions.mix(FINDPRIMERS.out.versions.first())

    CLUSTER ( FINDPRIMERS.out.filtered_fasta )
    ch_versions = ch_versions.mix(CLUSTER.out.versions.first())

    ch_clusters_ref = CLUSTER.out.consensus_fasta.combine(ch_ref, by:0)

    ALIGNCLUSTERS ( ch_clusters_ref )
    ch_versions = ch_versions.mix(ALIGNCLUSTERS.out.versions.first())

    SPLITBAM ( ALIGNCLUSTERS.out.bam )
    ch_versions = ch_versions.mix(SPLITBAM.out.versions.first())

    emit:
    fasta      = FINDPRIMERS.out.filtered_fasta  // channel: [ val(meta), [ fasta ] ]
    versions = ch_versions                       // channel: [ versions.yml ]
}

