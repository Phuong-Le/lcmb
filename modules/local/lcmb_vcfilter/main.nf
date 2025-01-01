process lcmbVcfilter {
    publishDir "${params.outdir}/filter_${mut_type}_out/${meta.pdid}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf)
    path(vcfilter_config)
    val mut_type

    output:
    tuple val(meta), path(vcf_filtered_gz), path(vcf_filtered_tbi)

    script:
    vcf_filtered = "${vcf.getName().tokenize(".").init().join(".")}.filter.vcf"
    vcf_filtered_gz = "${vcf_filtered}.gz"
    vcf_filtered_tbi = "${vcf_filtered_gz}.tbi"
    """
    vcfilter filter -o . -i ${vcfilter_config} ${vcf}
    bgzip ${vcf_filtered}
    tabix ${vcf_filtered_gz}
    """

}
