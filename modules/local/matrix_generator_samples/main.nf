process matrixGeneratorSamples {
    publishDir "${params.outdir}/filter_${mut_type}_out", mode: params.publish_dir_mode

    input:
    path vcf_ls, stageAs: 'sample_mutmat/*'
    val mut_type
    val sigprofiler_genome

    output:
    path mutmat_dir

    script:
    mutmat_dir = "sample_mutmat/output"
    """
    SigProfilerMatrixGenerator matrix_generator sample_mutmat ${sigprofiler_genome} sample_mutmat
    """

}
