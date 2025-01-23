process cgpVaf {
    label 'process_single_long'

    publishDir "${params.outdir}/filter_${mut_type}_out/${pdid}", mode: params.publish_dir_mode

    input:
    tuple val(pdid), val(sample_id_ls), val(match_normal_id), path(vcf), path(vcf_tbi), path(bam), path(bai), path(bas), path(met), path(bam_match), path(bai_match)
    val mut_type
    path fasta
    path fai
    path high_depth_regions
    path high_depth_regions_tbi

    output:
    tuple val(pdid), val(sample_id_ls), val(match_normal_id), path("${pdid}_${mut_type}_vaf.tsv")

    script:
    // meta = ['pdid':pdid, 'sample_id_ls':sample_id_ls, 'match_normal_id': match_normal_id]
    """
    rm -f "${pdid}_${mut_type}_vaf.tsv"
    rm -rf tmpvaf_*
    cgpVaf.pl -d . -o . -a ${mut_type} -g  ${fasta} -hdr  ${high_depth_regions} --vcf  ${vcf} --normal_bam ${bam_match} --tumour_bam ${bam} --normal_name ${match_normal_id} --tumour_name ${sample_id_ls.join(" ")} -ct 1

    grep -vE '^(##)' "${match_normal_id}_${sample_id_ls[0]}_${mut_type}_vaf.tsv" | sed 's/^#//' > "${pdid}_${mut_type}_vaf.tsv"
    """
}
