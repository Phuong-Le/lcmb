process conpairConcordanceUnmatchFilter {
    label 'process_tiny'

    publishDir "${params.outdir}/conpair_out", mode: params.publish_dir_mode

    input:
    path concordance
    path samples_path
    val concordance_threshold


    output:
    path input_filtered, emit: input_filtered
    path "*.log"

    script:
    input_filtered = "sample_paths_concordance_filtered.${samples_path.extension}"
    """
    conpair_concordance_filter_unmatch.py --samples_path ${samples_path} --concordance_path ${concordance} --concordance_threshold ${concordance_threshold} --outfile ${input_filtered}
    """
}
