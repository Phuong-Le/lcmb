process matrixGeneratorBranches {
    label 'process_tiny'

    publishDir "${params.outdir}/${outdir_basename}/${pdid}", mode: params.publish_dir_mode

    input:
    tuple val(pdid), path(muts_on_tree)
    val outdir_basename
    val sigprofiler_genome

    output:
    path matrix_generator_dir

    script:
    outdir = 'matrix_by_branch'
    matrix_generator_dir = "${outdir}/matrix_generator"
    """
    rm -rf ${outdir}
    mkdir ${outdir}
    split_vcf_to_branch.py --vcf_path ${muts_on_tree} --outdir ${outdir} --prefix ${pdid}
    SigProfilerMatrixGenerator matrix_generator ${pdid} ${sigprofiler_genome} ${matrix_generator_dir}
    rm ${matrix_generator_dir}/*.vcf
    """
}
