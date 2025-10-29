process rmPolyclonal {
    label 'process_tiny'

    publishDir "${params.outdir}/${outdir_basename}/${pdid}",  mode: params.publish_dir_mode


    input:
    tuple val(pdid), path(nr_path), path(nv_path), path(genotype_bin_path), path(clonality)
    val outdir_basename

    output:
    tuple val(pdid), path('NR_bbinom_filtered_clonal.txt'), path('NV_bbinom_filtered_clonal.txt'), path('genotype_bin_clonal.txt'), emit: phylogenetic_input


    script:
    """

    Rscript --vanilla ${projectDir}/bin/rm_polyclonal.R \
        --nr_path ${nr_path} \
        --nv_path ${nv_path} \
        --genotype_bin_path ${genotype_bin_path} \
        --clonality_path ${clonality} \
        --outdir .
    """
}
