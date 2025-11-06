process lcmbVcfilter {
    label 'process_tiny'

    publishDir "${params.outdir}/filter_${mut_type}_out/${meta.pdid}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf)
    path(vcfilter_config)
    val mut_type

    output:
    tuple val(meta), path(vcf_filtered_gz), path(vcf_filtered_tbi)

    script:
    // println(" lcmbFilter ${mut_type}: ${meta}   ${vcf}")
    vcf_filtered_gz = "${vcf.getName().tokenize(".").init().init().join(".")}.filter.vcf.gz"
    // vcf_filtered_gz = "${vcf_filtered}.gz"
    vcf_filtered_tbi = "${vcf_filtered_gz}.tbi"
    """
    vcfilter.sh --vcf ${vcf} --vcfilter_config ${vcfilter_config} --out_vcf_filtered ${vcf_filtered_gz}
    """

}
