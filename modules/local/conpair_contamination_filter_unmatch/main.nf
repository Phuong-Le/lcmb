process conpairContaminationUnmatchFilter {
    label 'process_tiny'

    publishDir "${params.outdir}/conpair_out", mode: params.publish_dir_mode

    input:
    path contamination
    path samples_path
    val contamination_threshold_samples


    output:
    path input_filtered, emit: input_filtered
    path "*.log"

    script:
    input_filtered = "sample_paths_concordance_filtered_contamination_filtered.${samples_path.extension}"
    """
    conpair_contamination_filter_unmatch.py --samples_path ${samples_path} --contamination_path ${contamination} --contamination_threshold_samples ${contamination_threshold_samples} --outfile ${input_filtered}
    """
}
