include { BOUNDARIES       } from '../../modules/local/boundaries'
include { FINDPRIMERS      } from '../../modules/local/findprimers'
include { CLUSTER          } from '../../modules/local/cluster'
include { ALIGNCLUSTERS    } from '../../modules/local/alignclusters'
include { TAGBAM           } from '../../modules/local/tagbam'
include { SPLITBAM         } from '../../modules/local/splitbam'
include { CALLCONSENSUS    } from '../../modules/local/callconsensus'
include { CIRCOS           } from '../../modules/local/circos'

workflow FILTERCLUSTERCONSENSUS {

    take:
    ch_samplesheet 
    ch_aligned_fasta   // channel: [ val(meta), [ fasta ] ]

    main:

    ch_primer1 = ch_samplesheet.map { meta, fastq, fasta, primer1, primer2, bed -> [meta, primer1] }
    ch_primer2 = ch_samplesheet.map { meta, fastq, fasta, primer1, primer2, bed -> [meta, primer2] }
    ch_ref = ch_samplesheet.map { meta, fastq, fasta, primer1, primer2, bed -> [meta, fasta] }
    ch_bed = ch_samplesheet.map { meta, fastq, fasta, primer1, primer2, bed -> [meta, bed] }
    ch_versions = Channel.empty()

    ch_fasta_primer1 = ch_aligned_fasta.combine(ch_primer1, by:0)
    ch_fasta_primer1_primer2 = ch_fasta_primer1.combine(ch_primer2, by:0)
    
    ch_ref_primer1 = ch_ref.combine(ch_primer1, by:0)
    ch_ref_primers = ch_ref_primer1.combine(ch_primer2, by:0)

    BOUNDARIES ( ch_ref_primers )
    ch_versions = ch_versions.mix(BOUNDARIES.out.versions.first())

    FINDPRIMERS ( ch_fasta_primer1_primer2 )
    ch_versions = ch_versions.mix(FINDPRIMERS.out.versions.first())

    CLUSTER ( FINDPRIMERS.out.filtered_fasta )
    ch_versions = ch_versions.mix(CLUSTER.out.versions.first())

    ch_clusters_ref = CLUSTER.out.consensus_fasta.combine(ch_ref, by:0)

    ALIGNCLUSTERS ( ch_clusters_ref )
    ch_versions = ch_versions.mix(ALIGNCLUSTERS.out.versions.first())

    TAGBAM ( ALIGNCLUSTERS.out.bam )
    ch_versions = ch_versions.mix(TAGBAM.out.versions.first())

    SPLITBAM ( TAGBAM.out.modified_bam )
    ch_versions = ch_versions.mix(SPLITBAM.out.versions.first())

    ch_bams_ref = SPLITBAM.out.bams.combine(ch_ref, by:0).transpose()
    
    CALLCONSENSUS ( ch_bams_ref )
    ch_versions = ch_versions.mix(CALLCONSENSUS.out.versions.first())

    // fastas = CALLCONSENSUS.out.con_fasta.map { file -> 
    //                 def key = file.name.toString().split('/').last().split('_clu').first()
    //                 return tuple(key, file) }.groupTuple()
    ch_vcfs = CALLCONSENSUS.out.vcf.map { file -> 
                    def key = file.name.toString().split('/').last().split('_clu').first()
                    return tuple(key, file) }.groupTuple()
    // ch_fastas_vcfs = fastas.combine(vcfs, by:0)
    //ch_fastas_vcfs.view()
    //ch_vcfs = CALLCONSENSUS.out.vcf
    ch_vcfs.view()
    ch_total_reads = SPLITBAM.out.txt
    ch_bounds = BOUNDARIES.out.txt
    ch_vcfs_bed = ch_vcfs.combine(ch_bed, by:0)
    ch_vcfs_bed_bounds = ch_vcfs_bed.combine(ch_bounds, by:0)
    ch_vcfs_bed_bounds_reads = ch_vcfs_bed_bounds.combine(ch_total_reads, by:0)
    ch_vcfs_bed_bounds_reads.view()
    // CIRCOS ( ch_vcfs_bed_bounds_reads )
    // ch_versions = ch_versions.mix(CIRCOS.out.versions.first())

    emit:
    fasta      = FINDPRIMERS.out.filtered_fasta  // channel: [ val(meta), [ fasta ] ]
    versions = ch_versions                       // channel: [ versions.yml ]
}

