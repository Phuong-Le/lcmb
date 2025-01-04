process betaBinomFilter {
    label 'process_tiny'

    publishDir "${params.outdir}/filter_${mut_type}_out/${pdid}", mode: params.publish_dir_mode

    input:
    tuple val(pdid), val(sample_id), val(match_normal_id), path(vcf), path(vcf_tbi), path(bed)
    val mut_type

    output:
    path vcf_filtered

    script:
    vcf_filtered = "${vcf.getName().tokenize(".").init().init().join(".")}.bbinom.vcf"
    """
    tabix -h -R ${bed} ${vcf} > ${vcf_filtered}
    """
}
