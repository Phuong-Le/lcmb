process betaBinomFilter {
    label 'process_tiny'

    publishDir "${params.outdir}/filter_${mut_type}_out/${pdid}", mode: params.publish_dir_mode

    input:
    tuple val(pdid), val(sample_id), val(match_normal_id), path(vcf), path(vcf_tbi), path(index_vcf)
    val mut_type
    path fai

    output:
    tuple val(pdid), val(sample_id), path(vcf_filtered), emit: vcf_filtered
    path vcf_filtered_gz, emit: vcf_filtered_gz
    path vcf_filtered_tbi, emit: vcf_filtered_tbi

    script:
    vcf_filtered = "${vcf.getName().tokenize(".").init().init().join(".")}.bbinom.vcf"
    vcf_filtered_gz = "${vcf_filtered}.gz"
    vcf_filtered_tbi = "${vcf_filtered_gz}.tbi"
    index_vcf_gz = "${index_vcf}.gz"
    """
    echo bgzipping
    bcftools view -O z -o ${index_vcf}_no_header.gz ${index_vcf}
    echo add header
    bcftools reheader -f ${fai} -o ${index_vcf}_with_header.gz ${index_vcf}_no_header.gz
    echo sorting
    bcftools sort -O z -o ${index_vcf_gz} ${index_vcf}_with_header.gz
    echo indexing
    bcftools index --tbi ${index_vcf_gz}

    bcftools isec -p info ${vcf} ${index_vcf_gz}

    cp info/0002.vcf ${vcf_filtered}
    bcftools view -O z -o ${vcf_filtered_gz} ${vcf_filtered}
    bcftools index --tbi ${vcf_filtered_gz}
    cp info/0002.vcf ${vcf_filtered}
    """
}
