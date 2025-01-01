/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { CONPAIR_FILTER_WITH_MATCH_NORMAL } from "${projectDir}/subworkflows/local/conpair_match"
include { LCMB_FILTER_SNV_MATCH } from "${projectDir}/subworkflows/local/lcmb_filter-snv_match"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LCMB {

    take:
    run_conpair
    run_filter_snv
    run_filter_indel
    run_phylogenetic
    ch_samplesheet // channel: samplesheet read in from --input
    input
    fasta
    fai
    fasta_dict
    marker_bed
    marker_txt
    concordance_threshold
    contamination_threshold_samples
    contamination_threshold_match
    snv_vcfilter_config
    snv_rho_threshold
    indel_vcfilter_config
    indel_rho_threshold
    high_depth_regions
    high_depth_regions_tbi
    hairpin_genome
    sigprofiler_genome

    main:

    if ( run_conpair == true ) {
        CONPAIR_FILTER_WITH_MATCH_NORMAL(
            ch_samplesheet,
            input,
            fasta,
            fai,
            fasta_dict,
            marker_bed,
            marker_txt,
            concordance_threshold,
            contamination_threshold_samples,
            contamination_threshold_match
        )
    }

    if ( run_filter_snv == true ) {
        if ( run_conpair = true ) {
            LCMB_FILTER_SNV_MATCH(
                CONPAIR_FILTER_WITH_MATCH_NORMAL.out
                .flatMap {
                samplesheetToList( it, "${projectDir}/assets/schemas/schema_input_conpair_filter-snv.json" )
                },
                snv_vcfilter_config,
                snv_rho_threshold,
                hairpin_genome,
                fasta,
                fai,
                high_depth_regions,
                high_depth_regions_tbi,
                sigprofiler_genome
            )
        }
        else {
            LCMB_FILTER_SNV_MATCH(
                ch_samplesheet,
                snv_vcfilter_config,
                snv_rho_threshold,
                hairpin_genome,
                fasta,
                fai,
                high_depth_regions,
                high_depth_regions_tbi,
                sigprofiler_genome
            )
        }
    }



}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
