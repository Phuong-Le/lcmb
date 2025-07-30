include { conpairPileup } from "$projectDir/modules/local/conpair_pileup"
include { verifyConcordance } from "$projectDir/modules/local/conpair_concordance"
include { conpairConcordanceUnmatchFilter } from "$projectDir/modules/local/conpair_concordance_filter_unmatch"
include { conpairContamination } from "$projectDir/modules/local/conpair_contamination"
include { conpairContaminationUnmatchFilter } from "$projectDir/modules/local/conpair_contamination_filter_unmatch"

workflow CONPAIR_FILTER_WITHOUT_MATCH_NORMAL {
    take:
    ch_samplesheet
    input
    fasta
    fai
    fasta_dict
    marker_bed
    marker_txt
    concordance_threshold
    contamination_threshold_samples
    contamination_threshold_match

    main:

    // pileup
    // sample
    conpairPileup(
        ch_samplesheet.
        map {
            meta, bam, bai, bam_match, bai_match ->
            tuple(meta.match_normal_id, meta.sample_id, bam, bai)
        },
        marker_bed,
        fasta,
        fai,
        fasta_dict
    )

    // Concordance
    pileup_pairs = conpairPileup.out
    .combine(conpairPileup.out)
    .filter {
        it[1] != it[4]
    }
    .map {
        match_normal_id1, sample_id1, pileup1, match_normal_id2, sample_id2, pileup2
        -> tuple(sample_id1, pileup1, sample_id2, pileup2)
    }

    concordance_all = verifyConcordance(
        pileup_pairs,
        marker_txt)
        .collectFile( name: 'concordance.txt', newLine: true, storeDir: "${params.outdir}/conpair_out" )


    // Filter by concordance
    conpairConcordanceUnmatchFilter(
        concordance_all,
        input,
        concordance_threshold
    )


    concordance_filtered_ch = conpairConcordanceUnmatchFilter.out.input_filtered
    .splitCsv(header: true, sep: '\t')
    .map { it -> tuple(it.sample_id, it.pdid) }
    .combine (
        conpairPileup.out
        .map {
            match_normal_id, sample_id, pileup
            -> tuple(sample_id, pileup)
        },
        by: 0
    )
    // .view{ "concordance_filtered_ch: $it" }

    concordance_filtered_pairs = concordance_filtered_ch
    .combine(concordance_filtered_ch, by: 1)
    .filter {
        it[1] != it[3]
    }
    .groupTuple(by: [0,1])
    .map {
        pdid, sample_id, pileups1, sample_ids2, pileups2
        -> tuple(sample_id, pileups1[0], sample_ids2[0], pileups2)
    }
    // .view { "contamination input: $it" }

    // Contamination
    contamination_all = conpairContamination(
        concordance_filtered_pairs,
        marker_txt)
        .collectFile( name: 'contamination.txt', newLine: true, storeDir: "${params.outdir}/conpair_out" )

    // Filtering contamination based on concordance and contamination
    conpairContaminationUnmatchFilter(
        contamination_all,
        input,
        contamination_threshold_samples
        )


    emit:
    conpairContaminationUnmatchFilter.out.input_filtered
}
