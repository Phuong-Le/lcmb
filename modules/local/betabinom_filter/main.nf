process betaBinomFilter {
    label 'process_tiny'

    publishDir "${params.outdir}/filter_${mut_type}_out/${pdid}", mode: params.publish_dir_mode

    input:
    tuple val(pdid), val(sample_id), val(match_normal_id), path(vcf), path(vcf_tbi), path(bed)
    val mut_type

    output:
    tuple val(pdid), val(sample_id), path(vcf_filtered), emit: vcf_filtered
    path vcf_filtered_gz, emit: vcf_filtered_gz
    path vcf_filtered_tbi, emit: vcf_filtered_tbi

    script:
    vcf_filtered = "${vcf.getName().tokenize(".").init().init().join(".")}.bbinom.vcf"
    vcf_filtered_gz = "${vcf_filtered}.gz"
    vcf_filtered_tbi = "${vcf_filtered_gz}.tbi"
    """
    tabix -h -R ${bed} ${vcf} > ${vcf_filtered}
    bgzip -c ${vcf_filtered} > ${vcf_filtered_gz}
    tabix -p vcf ${vcf_filtered_gz}
    """
}
