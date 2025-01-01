include { conpairPileup as conpairPileupSample } from "$projectDir/modules/local/conpair_pileup"
include { conpairPileup as conpairPileupMatch } from "$projectDir/modules/local/conpair_pileup"
include { verifyConcordance } from "$projectDir/modules/local/conpair_concordance"
include { conpairContamination } from "$projectDir/modules/local/conpair_contamination"
include { conpairFilter } from "$projectDir/modules/local/conpair_filter"

workflow CONPAIR_FILTER_WITH_MATCH_NORMAL {
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

    // check reference files specific to conpair exists
    fasta_dict = file(fasta_dict, checkIfExists: true)
    marker_txt = file(marker_txt, checkIfExists: true)
    marker_bed = file(marker_bed, checkIfExists: true)

    // pileup
    // sample
    conpairPileupSample(
        ch_samplesheet.
        map {
            meta, bam, bai, bam_match, bai_match ->
            tuple(meta.match_normal_id, meta.sample_id, bam, bai)
        },
        marker_bed,
        fasta,
        fai,
        fasta_dict)
    // normal
    conpairPileupMatch(
        ch_samplesheet.
        map {
            meta, bam, bai, bam_match, bai_match ->
            tuple(meta.match_normal_id, meta.match_normal_id, bam_match, bai_match)
        }.
        unique(),
        marker_bed,
        fasta,
        fai,
        fasta_dict)

    // Concordance between sample and match normal
    concordance_all = verifyConcordance(
        conpairPileupSample.out
        .combine(conpairPileupMatch.out)
        .map { sample -> tuple(sample[1], sample[2], sample[3], sample[5]) },
        marker_txt)
        .collectFile( name: 'conpair_out/concordance.txt', newLine: true )

    // Contamination
    contamination_all = conpairContamination(
        conpairPileupMatch.out
        .cross(conpairPileupSample.out)
        .map { sample -> tuple(sample[0][0], sample[0][2], sample[1][1], sample[1][2]) },
        marker_txt)
        .collectFile( name: 'conpair_out/contamination.txt', newLine: true )

    // Filtering contamination based on concordance and contamination
    conpairFilter(
        concordance_all,
        contamination_all,
        input,
        concordance_threshold,
        contamination_threshold_samples,
        contamination_threshold_match)


    emit:
    conpairFilter.out.input_filtered
}
