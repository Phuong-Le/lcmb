#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/lcmb
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/lcmb
    Website: https://nf-co.re/lcmb
    Slack  : https://nfcore.slack.com/channels/lcmb
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { LCMB  } from './workflows/lcmb'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_lcmb_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_lcmb_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_lcmb_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.fasta = getGenomeAttribute('fasta')
params.fai = getGenomeAttribute('fai')
params.fasta_dict = getGenomeAttribute('fasta_dict')
params.marker_bed = getGenomeAttribute('marker_bed')
params.marker_txt = getGenomeAttribute('marker_txt')
params.high_depth_regions = getGenomeAttribute('high_depth_regions')
params.high_depth_regions_tbi = getGenomeAttribute('high_depth_regions_tbi')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_LCMB {

    take:
    run_conpair
    run_filter_snv
    run_filter_indel
    run_phylogenetics
    samplesheet // channel: samplesheet read in from --input
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

    //
    // WORKFLOW: Run pipeline
    //
    LCMB (
        run_conpair,
        run_filter_snv,
        run_filter_indel,
        run_phylogenetics,
        samplesheet,
        input,
        fasta,
        fai,
        fasta_dict,
        marker_bed,
        marker_txt,
        concordance_threshold,
        contamination_threshold_samples,
        contamination_threshold_match,
        snv_vcfilter_config,
        snv_rho_threshold,
        indel_vcfilter_config,
        indel_rho_threshold,
        high_depth_regions,
        high_depth_regions_tbi,
        hairpin_genome,
        sigprofiler_genome
    )


}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.run_conpair,
        params.run_filter_snv,
        params.run_filter_indel,
        params.run_phylogenetics
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_LCMB (
        params.run_conpair,
        params.run_filter_snv,
        params.run_filter_indel,
        params.run_phylogenetics,
        PIPELINE_INITIALISATION.out.samplesheet,
        params.input,
        params.fasta,
        params.fai,
        params.fasta_dict,
        params.marker_bed,
        params.marker_txt,
        params.concordance_threshold,
        params.contamination_threshold_samples,
        params.contamination_threshold_match,
        params.snv_vcfilter_config,
        params.snv_rho_threshold,
        params.indel_vcfilter_config,
        params.indel_rho_threshold,
        params.high_depth_regions,
        params.high_depth_regions_tbi,
        params.hairpin_genome,
        params.sigprofiler_genome
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        null
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
