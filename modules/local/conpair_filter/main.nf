process conpairFilter {
    label 'process_tiny'

    publishDir "${params.outdir}/conpair_out", mode: params.publish_dir_mode

    input:
    path concordance
    path contamination
    path samples_path
    val concordance_threshold
    val contamination_threshold_samples
    val contamination_threshold_match


    output:
    path input_filtered, emit: input_filtered
    path "*.log"

    script:
    input_filtered = "sample_paths_contamination_filtered.${samples_path.extension}"
    """
    conpair_contamination_filter.py --samples_path ${samples_path} --concordance_path ${concordance} --contamination_path ${contamination} --concordance_threshold ${concordance_threshold} --contamination_threshold_samples ${contamination_threshold_samples} --contamination_threshold_match ${contamination_threshold_match} --outfile ${input_filtered}
    """
}
