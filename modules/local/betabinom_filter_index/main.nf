process betaBinomFilterIndex {
    label 'process_tiny'

    publishDir "${params.outdir}/filter_${mut_type}_out/${pdid}", mode: params.publish_dir_mode

    // Beta Binomial filtering of germline mutations and artefacts, based on Tim Coorens' R script
    // The outcome is a bed file for the indices of the PASSED mutations
    input:
    tuple val(pdid), val(sample_id_ls), val(match_normal_id), path(vaf)
    val mut_type
    val rho_threshold

    output:
    tuple val(pdid), path("*.bed"), emit: betabinom_bed
    path "germline_ids.txt", optional: true, emit: germline_ids
    path "somatic_ids.txt", optional: true, emit: somatic_ids
    path "somatic_ids_rho.txt", optional: true, emit: somatic_ids_rho
    tuple val(pdid), path("NR_bbinom_filtered.txt"), path("NV_bbinom_filtered.txt"), path("genotype_bin.txt"), optional: true, emit: phylogenetic_raw_input

    script:
    """
    Rscript --vanilla ${projectDir}/bin/beta_binom_filter_index.R --libpath=${projectDir}/bin --cgpvaf_out=${vaf} --match_normal_id=${match_normal_id} --rho_threshold=${rho_threshold} --outdir=.
    """
}
