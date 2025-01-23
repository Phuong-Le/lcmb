process concatMatrices {
    label 'process_tiny'

    publishDir "${params.outdir}/${outdir_basename}", mode: params.publish_dir_mode

    input:
    path matrix_dirs, stageAs: "indir/matrix_generator*"
    val outdir_basename

    output:
    path outdir

    script:
    outdir = "combined_matrices_by_branches"
    """
    rm -rf ${outdir}
    mkdir ${outdir}
    concat_mutmats.py --indir indir --outdir ${outdir}
    """

}
