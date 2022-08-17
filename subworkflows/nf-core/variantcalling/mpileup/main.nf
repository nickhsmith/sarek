include { BCFTOOLS_MPILEUP                       } from '../../../../modules/nf-core/modules/bcftools/mpileup/main'
include { GATK4_MERGEVCFS   as MPILEUP_MERGEVCFS } from '../../../../modules/nf-core/modules/gatk4/mergevcfs/main'

workflow RUN_MPILEUP {
    take:
        cram                    // channel: [mandatory] [meta, cram, interval]
        dict                    // channel: [mandatory]
        fasta                   // channel: [mandatory]

    main:

    ch_versions = Channel.empty()

    BCFTOOLS_MPILEUP(cram, fasta, false)
    mpileup = BCFTOOLS_MPILEUP.out.vcf.branch{
            intervals:    it[0].num_intervals > 1  //check meta.num_intervals
            no_intervals: it[0].num_intervals <= 1 //check meta.num_intervals
        }

    //Merge mpileup only when intervals and natural order sort them
    MPILEUP_MERGEVCFS(mpileup.intervals
            .map{ meta, vcf ->
                new_meta = meta.tumor_id ? [
                                                id:             meta.tumor_id + "_vs_" + meta.normal_id,
                                                normal_id:      meta.normal_id,
                                                num_intervals:  meta.num_intervals,
                                                patient:        meta.patient,
                                                sex:            meta.sex,
                                                tumor_id:       meta.tumor_id,
                                            ] // not annotated, so no variantcaller necessary
                                            : [
                                                id:             meta.sample,
                                                num_intervals:  meta.num_intervals,
                                                patient:        meta.patient,
                                                sample:         meta.sample,
                                                status:         meta.status,
                                                sex:            meta.sex,
                                            ]
                [groupKey(new_meta, meta.num_intervals), vcf]
            }.groupTuple(sort:true),
        dict
    )


    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)
    ch_versions = ch_versions.mix(MPILEUP_MERGEVCFS.out.versions)

    emit:
    versions = ch_versions
    mpileup = Channel.empty().mix(MPILEUP_MERGEVCFS.out.vcf, mpileup.no_intervals)
}
