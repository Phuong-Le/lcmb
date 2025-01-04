process verifyConcordance {
    label 'process_tiny'


    input:
    tuple val(sample_id), path(sample_pileup), val(match_normal_id), path(match_pileup)
    path(marker_txt)

    output:
    path(concordance_path)

    script:
    concordance_path = "${sample_id}_${match_normal_id}.concordance.txt"

    """
    verify_concordance.py -T ${sample_pileup} -N ${match_pileup} -M ${marker_txt} -O ${concordance_path}
    """
}
