process mutToTree {
    label 'process_tiny'

    publishDir "${params.outdir}/${outdir_basename}/${pdid}", mode: params.publish_dir_mode

    input:
    tuple val(pdid), path(nr_path), path(nv_path), path(genotype_bin_path), path(topology)
    val outdir_basename

    output:
    tuple val(pdid), path(muts_on_tree_txt), emit: muts_on_tree
    tuple val(pdid), path(tree_w_branchlength), path(tree_w_branchlength_plot)

    script:
    tree_w_branchlength = "${pdid}.tree_with_branch_length.tree"
    tree_w_branchlength_plot = "${pdid}.tree_with_branch_length.pdf"
    muts_on_tree_txt = "${pdid}.muts_assigned_to_tree.txt"
    """
    Rscript --vanilla ${projectDir}/bin/assign_mut_to_tree.R --nr_path=$nr_path --nv_path=$nv_path --genotype_bin_path=$genotype_bin_path --tree_path=$topology --out_prefix=$pdid
    """

}
