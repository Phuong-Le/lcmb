process getPhylogeny {
    label 'process_single'

    publishDir "${params.outdir}/${outdir_basename}/${pdid}", mode: params.publish_dir_mode

    input:
    tuple val(pdid), path(nr_path), path(nv_path), path(genotype_bin_path)
    val outdir_basename

    output:
    tuple val(pdid), path(muts_on_tree_txt), emit: muts_on_tree
    tuple val(pdid), path(tree_path), emit: topology
    tuple val(pdid), path(fasta), path(tree_w_branchlength), path(tree_w_branchlength_plot)
    path mpboot_log

    script:
    fasta = "${pdid}.fasta"
    tree_path = "${fasta}.treefile"
    mpboot_log = "${fasta}.log"
    tree_w_branchlength = "${pdid}.tree_with_branch_length.tree"
    tree_w_branchlength_plot = "${pdid}.tree_with_branch_length.pdf"
    muts_on_tree_txt = "${pdid}.muts_assigned_to_tree.txt"
    """
    genotype_bin_to_fasta.py --genotype_bin_path ${genotype_bin_path} --outfile ${fasta}
    mpboot -s ${fasta} -bb 1000 -cost e -st DNA
    Rscript --vanilla ${projectDir}/bin/assign_mut_to_tree.R --nr_path=$nr_path --nv_path=$nv_path --genotype_bin_path=$genotype_bin_path --tree_path=$tree_path --out_prefix=$pdid
    """

}
