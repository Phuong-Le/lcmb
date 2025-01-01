process spectraPlottingSamples {
    publishDir "${params.outdir}/filter_${mut_type}_out", mode: params.publish_dir_mode

    input:
    path mutmat_dir, stageAs: 'sample_mutmat/*'
    val mut_type

    output:
    path mutmat_plot_dir

    script:
    mutmat_plot_dir = "sample_mutmat/plots"
    """
    sigprofiler_plotting.py --matrix_dir ${mutmat_dir} --project sample_mutmat --output_path ${mutmat_plot_dir}
    """
}
