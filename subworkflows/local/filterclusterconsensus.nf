include { BOUNDARIES       } from '../../modules/local/boundaries'
include { FINDPRIMERS      } from '../../modules/local/findprimers'
include { CLUSTER          } from '../../modules/local/cluster'
include { ALIGNCLUSTERS    } from '../../modules/local/alignclusters'
include { TAGBAM           } from '../../modules/local/tagbam'
include { SPLITBAM         } from '../../modules/local/splitbam'
include { CALLCONSENSUS    } from '../../modules/local/callconsensus'
include { CALLCONSENSUSPP  } from '../../modules/local/callconsensuspp'
include { SPECIESPLOTS     } from '../../modules/local/speciesplots/main'
include { SPECIESPLOTNOPP } from '../../modules/local/speciesplotnopp/main'
workflow FILTERCLUSTERCONSENSUS {

    take:
    ch_samplesheet 
    ch_aligned_fasta   // channel: [ val(meta), [ fasta ] ]

    main:

    ch_primer1 = ch_samplesheet.map { meta, fastq, fasta, primer1, primer2, bed, introns -> [meta, primer1] }
    ch_primer2 = ch_samplesheet.map { meta, fastq, fasta, primer1, primer2, bed, introns -> [meta, primer2] }
    ch_ref = ch_samplesheet.map { meta, fastq, fasta, primer1, primer2, bed, introns -> [meta, fasta] }
    ch_bed = ch_samplesheet.map { meta, fastq, fasta, primer1, primer2, bed, introns -> 
                                                                            meta = meta.id.split('_').last()
                                                                            def key = meta
                                                                            return tuple(key, bed) }.groupTuple().map {group -> 
                                                                                                                def (key, values) = group
                                                                                                                [key, values[0]]} 
    ch_fastq = ch_samplesheet.map { meta, fastq, fasta, primer1, primer2, bed, introns -> [meta, fastq] }
    ch_extra_fasta = ch_fastq.map { meta, file -> 
                    def fileType = file.name.toString().split('/').last().split('\\.').last()
                    if (fileType == "fasta") {
                        return tuple(meta, file)
                    } 
                 }
    ch_extra_fasta_ref = ch_extra_fasta.combine(ch_ref, by:0)

    ch_versions = Channel.empty()

    ch_fasta_primer1 = ch_aligned_fasta.combine(ch_primer1, by:0)
    ch_fasta_primer1_primer2 = ch_fasta_primer1.combine(ch_primer2, by:0)
    
    ch_ref_primer1 = ch_ref.combine(ch_primer1, by:0)
    ch_ref_primers = ch_ref_primer1.combine(ch_primer2, by:0)

    BOUNDARIES ( ch_ref_primers )
    ch_versions = ch_versions.mix(BOUNDARIES.out.versions.first())

    FINDPRIMERS ( ch_fasta_primer1_primer2 )
    ch_versions = ch_versions.mix(FINDPRIMERS.out.versions.first())

    ch_find_primers_bounds = FINDPRIMERS.out.filtered_fasta.combine(BOUNDARIES.out.txt, by:0)
    ch_find_primers_bounds.view()

    CLUSTER ( ch_find_primers_bounds )
    ch_versions = ch_versions.mix(CLUSTER.out.versions.first())

    ch_clusters_ref = CLUSTER.out.consensus_fasta.combine(ch_ref, by:0)
    // Here combine the custom user added fasta file (of known preprocessed pseudogenes in some cases) with the other samples if it exixts

    ch_clusters_ref = ch_clusters_ref.concat(ch_extra_fasta_ref)

    ALIGNCLUSTERS ( ch_clusters_ref )
    ch_versions = ch_versions.mix(ALIGNCLUSTERS.out.versions.first())

    TAGBAM ( ALIGNCLUSTERS.out.bam )
    ch_versions = ch_versions.mix(TAGBAM.out.versions.first())

    SPLITBAM ( TAGBAM.out.modified_bam )
    ch_versions = ch_versions.mix(SPLITBAM.out.versions.first())

    ch_bams_ref = SPLITBAM.out.bams.combine(ch_ref, by:0).transpose()
    // Need to split out the processed pseudogene bams to go into separate call consensus modules, then merge for input into circos module
    ch_bams_ref_pp = ch_bams_ref.map { meta, file, ref -> 
                    def fileType = file.name.toString().split('/').last().split('\\.bam').first().split('_').last()
                    if (fileType.contains("pp")) {
                        return tuple(meta, file, ref)
                    } 
                 }
     ch_bams_ref = ch_bams_ref.map { meta, file, ref -> 
                    def fileType = file.name.toString().split('/').last().split('\\.bam').first().split('_').last()
                    if (!fileType.contains("pp")) {
                        return tuple(meta, file, ref)
                    } 
                 }
    CALLCONSENSUS ( ch_bams_ref )
    ch_versions = ch_versions.mix(CALLCONSENSUS.out.versions.first())

    CALLCONSENSUSPP ( ch_bams_ref_pp )
    ch_versions = ch_versions.mix(CALLCONSENSUSPP.out.versions.first())

    // combine all samples that queried the same gene so can compare across samples
    ch_orf_beds_not_ref = CALLCONSENSUS.out.orf_bed_tr.map { file -> 
                    def key = file.name.toString().split('/').last().split('_clu').first().split('_').last()
                    return tuple(key, file) }.groupTuple()

    ch_orf_beds_pp_not_ref = CALLCONSENSUSPP.out.orf_bed_tr.map { file -> 
                    def key = file.name.toString().split('/').last().split('_pp').first().split('_').last()
                    return tuple(key, file) }.groupTuple()

    ch_vcfs = CALLCONSENSUS.out.vcf.map { file -> 
                    def key = file.name.toString().split('/').last().split('_clu').first().split('_').last()
                    return tuple(key, file) }.groupTuple()

    ch_vcfs_pp = CALLCONSENSUSPP.out.vcf.map { file -> 
                def key = file.name.toString().split('/').last().split('_pp').first().split('_').last()
                return tuple(key, file) }.groupTuple()

    ch_total_reads = SPLITBAM.out.txt.map { meta, txt -> 
                                    meta = meta.id.split('_').last()
                                    def key = meta
                                    return tuple(key, txt) }.groupTuple()

    ch_vcfs_all = ch_vcfs.combine(ch_vcfs_pp, by:0).map {meta, clusters, pps -> 
                                                            def files = clusters + pps
                                                            return tuple(meta, files) }

   
    ch_orf_beds_all = ch_orf_beds_not_ref.combine(ch_orf_beds_pp_not_ref, by:0).map {meta, clusters, pps ->
                                                            def files = clusters + pps
                                                            return tuple(meta, files)}

    ch_bounds = BOUNDARIES.out.txt.map { meta, txt -> 
                                    meta = meta.id.split('_').last()
                                    def key = meta
                                    return tuple(key, txt) }.groupTuple().map {group -> 
                                                                        def (key, values) = group
                                                                        [key, values[0]]}

    ch_vcfs_bed = ch_vcfs_all.combine(ch_bed, by:0)
    ch_vcfs_bed_bounds = ch_vcfs_bed.combine(ch_bounds, by:0)
    ch_vcfs_bed_bounds_reads = ch_vcfs_bed_bounds.combine(ch_total_reads, by:0)
    ch_vcfs_bed_bounds_reads_orfs = ch_vcfs_bed_bounds_reads.combine(ch_orf_beds_all, by:0)
    
    SPECIESPLOTS ( ch_vcfs_bed_bounds_reads_orfs )
    ch_versions = ch_versions.mix(SPECIESPLOTS.out.versions.first())

    // if no processed pseudogenes, then run separately 
    ch_vcfs_bed2 = ch_vcfs.combine(ch_bed, by:0)
    ch_vcfs_bed_bounds2 = ch_vcfs_bed2.combine(ch_bounds, by:0)
    ch_vcfs_bed_bounds_reads2 = ch_vcfs_bed_bounds2.combine(ch_total_reads, by:0)
    ch_vcfs_bed_bounds_reads_orfs2 = ch_vcfs_bed_bounds_reads2.combine(ch_orf_beds_not_ref, by:0)

    SPECIESPLOTNOPP ( ch_vcfs_bed_bounds_reads_orfs2 )
    ch_versions = ch_versions.mix(SPECIESPLOTNOPP.out.versions.first())

    emit:
    fasta    = FINDPRIMERS.out.filtered_fasta    // channel: [ val(meta), [ fasta ] ]
    versions = ch_versions                       // channel: [ versions.yml ]
}

