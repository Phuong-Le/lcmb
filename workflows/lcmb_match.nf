/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { CONPAIR_FILTER_WITH_MATCH_NORMAL } from "${projectDir}/subworkflows/local/conpair_match"
include { LCMB_FILTER_SNV_MATCH } from "${projectDir}/subworkflows/local/lcmb_filter-snv_match"
include { LCMB_FILTER_INDEL_MATCH } from "${projectDir}/subworkflows/local/lcmb_filter-indel_match"
include { PHYLOGENETICS } from "${projectDir}/subworkflows/local/phylogenetics"
include { PHYLOGENETICS_PROVIDED_TREE_TOPOLOGY } from "${projectDir}/subworkflows/local/phylogenetics_provided_topology"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LCMB_MATCH {

    take:
    run_conpair
    run_filter_snv
    run_filter_indel
    run_phylogenetics
    ch_samplesheet_conpair // channel: samplesheet read in from --input
    ch_samplesheet_filter_snv
    ch_samplesheet_filter_indel
    ch_samplesheet_phylogenetics
    ch_samplesheet_snv_then_indel
    ch_samplesheet_topology
    ch_samplesheet_clonality
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
    hairpin2_input_snv_json
    hairpin2_name_mapping_snv
    hairpin2_input_indel_json
    hairpin2_name_mapping_indel
    sigprofiler_genome
    min_good_reads
    max_K
    max_iter
    nchains
    clonal_threshold
    proportion_pass_clonality
    snv_then_indel
    provided_topology
    phylogenetics_outdir_basename

    main:

    // CONPAIR
    if ( run_conpair == true ) {
        CONPAIR_FILTER_WITH_MATCH_NORMAL(
            ch_samplesheet_conpair,
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

    // FILTER SNV
    if ( run_filter_snv == true ) {
        if ( run_conpair == true ) {
            // LCMB_FILTER_SNV_MATCH(
            //     CONPAIR_FILTER_WITH_MATCH_NORMAL.out
            //     .flatMap {
            //     samplesheetToList( it, "${projectDir}/assets/schemas/schema_input_conpair_filter-snv.json" )
            //     },
            //     snv_vcfilter_config,
            //     snv_rho_threshold,
            //     hairpin_genome,
            //     fasta,
            //     fai,
            //     high_depth_regions,
            //     high_depth_regions_tbi,
            //     sigprofiler_genome
            // )

            LCMB_FILTER_SNV_MATCH(
                CONPAIR_FILTER_WITH_MATCH_NORMAL.out
                .splitCsv( header: true, sep : '\t' )
                .map { row -> tuple( ['sample_id': row.sample_id, 'match_normal_id': row.match_normal_id, 'pdid' : row.pdid], row.bam, row.bai, row.bas, row.met, row.bam_match, row.bai_match, row.snv_vcf, row.snv_vcf_tbi ) },
                snv_vcfilter_config,
                snv_rho_threshold,
                hairpin2_input_snv_json,
                hairpin2_name_mapping_snv,
                fasta,
                fai,
                high_depth_regions,
                high_depth_regions_tbi,
                sigprofiler_genome,
                min_good_reads,
                max_K,
                max_iter,
                nchains,
                clonal_threshold,
                proportion_pass_clonality
            )


        }
        else {
            LCMB_FILTER_SNV_MATCH(
                ch_samplesheet_filter_snv,
                snv_vcfilter_config,
                snv_rho_threshold,
                hairpin2_input_snv_json,
                hairpin2_name_mapping_snv,
                fasta,
                fai,
                high_depth_regions,
                high_depth_regions_tbi,
                sigprofiler_genome,
                min_good_reads,
                max_K,
                max_iter,
                nchains,
                clonal_threshold,
                proportion_pass_clonality
            )
            // LCMB_FILTER_SNV_MATCH.out.view()
        }
    }

    // FILTER INDEL
    if ( run_filter_indel == true ) {
        if ( run_conpair == true ) {
            // LCMB_FILTER_INDEL_MATCH(
            //     CONPAIR_FILTER_WITH_MATCH_NORMAL.out
            //     .flatMap {
            //     samplesheetToList( it, "${projectDir}/assets/schemas/schema_input_conpair_filter-indel.json" )
            //     },
            //     indel_vcfilter_config,
            //     indel_rho_threshold,
            //     fasta,
            //     fai,
            //     high_depth_regions,
            //     high_depth_regions_tbi,
            //     sigprofiler_genome
            // )


            LCMB_FILTER_INDEL_MATCH(
                CONPAIR_FILTER_WITH_MATCH_NORMAL.out
                .splitCsv( header: true, sep : '\t' )
                .map { row -> tuple( ['sample_id': row.sample_id, 'match_normal_id': row.match_normal_id, 'pdid' : row.pdid], row.bam, row.bai, row.bas, row.met, row.bam_match, row.bai_match, row.indel_vcf, row.indel_vcf_tbi ) },
                indel_vcfilter_config,
                indel_rho_threshold,
                hairpin2_input_indel_json,
                hairpin2_name_mapping_indel,
                fasta,
                fai,
                high_depth_regions,
                high_depth_regions_tbi,
                sigprofiler_genome
            )
        }
        else {
            LCMB_FILTER_INDEL_MATCH(
                ch_samplesheet_filter_indel,
                indel_vcfilter_config,
                indel_rho_threshold,
                hairpin2_input_indel_json,
                hairpin2_name_mapping_indel,
                fasta,
                fai,
                high_depth_regions,
                high_depth_regions_tbi,
                sigprofiler_genome
            )
        }
    }

    // PHYLOGENETICS
    if ( run_phylogenetics == true ) {
        if ( run_filter_snv == true ) {
            // only run this if there are more than 2 sample per donor (genotype_bin only has one column)
            // phylogenetics without tree topology provided
            PHYLOGENETICS(
                LCMB_FILTER_SNV_MATCH.out
                .filter { it[3].readLines().first().split(' ').size() > 2 },
                'phylogenetics_snp_out',
                sigprofiler_genome
            )
            if ( run_filter_indel == true ) {
                // only run this if there are more than 2 sample per donor (genotype_bin only has one column)

                topology_ch = PHYLOGENETICS.out
                    .filter { it != null }
                // get clonality from SNV filtering
                clonality_ch = LCMB_FILTER_SNV_MATCH.out
                    .filter { it[3].readLines().first().split(' ').size() > 2 }
                    .map { pdid, nr_path, nv_path, genotype_bin_path, clonality ->
                        tuple(pdid, clonality)
                    }

                // phylogenetics provided tree topology
                PHYLOGENETICS_PROVIDED_TREE_TOPOLOGY(
                    LCMB_FILTER_INDEL_MATCH.out
                    .filter { it[3].readLines().first().split(' ').size() > 2 }
                    .combine( clonality_ch, by: 0)
                    .combine( topology_ch, by: 0 ),
                    'phylogenetics_indel_out',
                    sigprofiler_genome
                )
            }
        }
        else if ( run_filter_indel == true ) {
            // only run this if there are more than 2 sample per donor (genotype_bin only has one column)
            // phylogenetics provided tree topology
            topology_ch = ch_samplesheet_topology
                .map {
                    meta, topology -> tuple(meta.pdid, topology)
                }
                .filter { it != null }
                .unique()
            clonality_ch = ch_samplesheet_clonality
                .map {
                    meta, clonality -> tuple(meta.pdid, clonality)
                }
                .filter { it != null }
                .unique()

            PHYLOGENETICS_PROVIDED_TREE_TOPOLOGY(
                LCMB_FILTER_INDEL_MATCH.out
                    .filter { it[3].readLines().first().split(' ').size() > 2 }
                    .combine(
                        clonality_ch,
                        by: 0
                    )
                    .combine(
                        topology_ch,
                        by: 0
                    ),
                'phylogenetics_indel_out',
                sigprofiler_genome
            )
        }
        else {
            // phylogenetics only

            clonality_ch = ch_samplesheet_clonality
                .map {
                    meta, clonality -> tuple(meta.pdid, clonality)
                }
                .filter { it != null }
                .unique()
            // pipeline selection
            if ( snv_then_indel == true ) {
                // phylogenetics for SNV
                PHYLOGENETICS(
                    ch_samplesheet_snv_then_indel
                    .map {
                        meta, nr_path_snv, nv_path_snv, genotype_bin_path_snv, nr_path_indel, nv_path_indel, genotype_bin_path_indel
                        -> tuple(meta.pdid, nr_path_snv, nv_path_snv, genotype_bin_path_snv)
                    }
                    .filter { it[3].readLines().first().split(' ').size() > 2 }
                    .combine( clonality_ch, by: 0 )
                    , 'phylogenetics_snp_out',
                    sigprofiler_genome
                )
                // phylogenetics for INDEL
                PHYLOGENETICS_PROVIDED_TREE_TOPOLOGY(
                    ch_samplesheet_snv_then_indel
                    .map {
                            meta, nr_path_snv, nv_path_snv, genotype_bin_path_snv, nr_path_indel, nv_path_indel, genotype_bin_path_indel
                            -> tuple(meta.pdid, nr_path_indel, nv_path_indel, genotype_bin_path_indel)
                        }
                    .filter { it[3].readLines().first().split(' ').size() > 2 }
                    .combine( clonality_ch, by: 0 )
                    .combine( PHYLOGENETICS.out, by: 0),
                    'phylogenetics_indel_out',
                    sigprofiler_genome
                )

            }
            else if ( provided_topology == true ) {
                phylogenetics_outdir_basename = (phylogenetics_outdir_basename == null) ? 'phylogenetics_indel_out' : phylogenetics_outdir_basename

                PHYLOGENETICS_PROVIDED_TREE_TOPOLOGY(
                    ch_samplesheet_phylogenetics
                    .map {
                        meta, nr_path, nv_path, genotype_bin_path
                        -> tuple(meta.pdid, nr_path, nv_path, genotype_bin_path)
                    }
                    .filter { it[3].readLines().first().split(' ').size() > 2 }
                    .combine( clonality_ch, by: 0 )
                    .combine( ch_samplesheet_topology, by: 0 ),
                    phylogenetics_outdir_basename,
                    sigprofiler_genome
                )
            }
            else {
                phylogenetics_outdir_basename = (phylogenetics_outdir_basename == null) ? 'phylogenetics_snp_out' : phylogenetics_outdir_basename

                PHYLOGENETICS(
                    ch_samplesheet_phylogenetics
                    .map {
                        meta, nr_path, nv_path, genotype_bin_path
                        -> tuple(meta.pdid, nr_path, nv_path, genotype_bin_path)
                    }
                    .filter { it[3].readLines().first().split(' ').size() > 2 }
                    .combine( clonality_ch, by: 0 ),
                    phylogenetics_outdir_basename,
                    sigprofiler_genome
                )
            }
        }
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
