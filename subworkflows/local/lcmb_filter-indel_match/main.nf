include { lcmbVcfilter } from "$projectDir/modules/local/lcmb_vcfilter"
include { getChromCgpVaf } from "$projectDir/modules/local/get_chrom_cgpvaf"
include { cgpVafChrom } from "$projectDir/modules/local/cgpvaf_chrom"
include { cgpVafConcat } from "$projectDir/modules/local/cgpvaf_concat"
include { betaBinomFilterIndex } from "$projectDir/modules/local/betabinom_filter_index"
include { betaBinomFilter } from "$projectDir/modules/local/betabinom_filter"
include { matrixGeneratorSamples } from "$projectDir/modules/local/matrix_generator_samples"
include { spectraPlottingSamples } from "$projectDir/modules/local/spectra_plotting_samples"



workflow LCMB_FILTER_INDEL_MATCH {
    take:
    input
    vcfilter_config
    rho_threshold
    fasta
    fai
    high_depth_regions
    high_depth_regions_tbi
    sigprofiler_genome

    main:

    // setup
    mut_type = 'indel'
    bams = input
        .map {
            meta, bam, bai, bas, met, bam_match, bai_match, vcf, vcf_tbi ->
            tuple(meta, bam, bai, bas, met, bam_match, bai_match)
        }

    // LCMB vcfilter
    lcmbVcfilter(
        input.
        map {
            meta, bam, bai, bas, met, bam_match, bai_match, vcf, vcf_tbi ->
            tuple(meta, vcf)
        },
        vcfilter_config,
        mut_type
    )

    // cgpVaf
    getChromCgpVaf(
        fai,
        high_depth_regions
    )

    vaf_input_files = lcmbVcfilter.out
        .combine( bams, by: 0 )
        .map {
            meta, vcf_filtered_gz, vcf_filtered_tbi, bam, bai, bas, met, bam_match, bai_match ->
            tuple( meta.pdid, meta.sample_id, meta.match_normal_id, vcf_filtered_gz, vcf_filtered_tbi, bam, bai, bas, met, bam_match, bai_match )
        }
        .groupTuple ( by: 0 )
        .map {
            pdid, sample_id, match_normal_id, vcf_filtered_gz, vcf_filtered_tbi, bam, bai, bas, met, bam_match, bai_match
            -> tuple(pdid, sample_id, match_normal_id[0], vcf_filtered_gz, vcf_filtered_tbi, bam, bai, bas, met, bam_match[0], bai_match[0])
        }

    cgpVafChrom(
        vaf_input_files
        .combine(
            getChromCgpVaf.out
            .splitCsv(header: false, sep: "\t")
        ),
        mut_type,
        fasta,
        fai,
        high_depth_regions,
        high_depth_regions_tbi
    )

    cgpVafConcat(
        cgpVafChrom.out
        .groupTuple( by: [0, 1, 2] )
        .combine(
            vaf_input_files,
            by: [0, 1, 2]
        ),
        mut_type,
        fasta,
        fai,
        high_depth_regions,
        high_depth_regions_tbi
    )

    // BetaBinomial filtering for germline and LCM artefacts based on cgpVaf (methods by Tim Coorens)
    betaBinomFilterIndex(
        cgpVafConcat.out,
        mut_type,
        rho_threshold
    )

    // betabinomial filtering on the vcfiltered vcfs
    betaBinomFilter(
        lcmbVcfilter.out
        .map {
            meta, vcf_filtered_gz, vcf_filtered_tbi ->
            tuple(meta.pdid, meta.sample_id, meta.match_normal_id, vcf_filtered_gz, vcf_filtered_tbi)
        }
        .combine(
            betaBinomFilterIndex.out.betabinom_bed,
            by: 0
        ),
        mut_type
    )

    // generate mutation matrix for the samples by SigProfilerMatrixGenerator
    matrixGeneratorSamples(
        betaBinomFilter.out.vcf_filtered.toList(),
        mut_type,
        sigprofiler_genome
        )

    // plot spectra
    spectraPlottingSamples(
        matrixGeneratorSamples.out,
        mut_type
        )

    emit:
    betaBinomFilterIndex.out.phylogenetics_input


}
