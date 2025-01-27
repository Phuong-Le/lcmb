process cgpVafChrom {
    label 'process_single'

    input:
    tuple val(pdid), val(sample_id_ls), val(match_normal_id), path(vcf), path(vcf_tbi), path(bam), path(bai), path(bas), path(met), path(bam_match), path(bai_match), val(chrom)
    val mut_type
    path fasta
    path fai
    path high_depth_regions
    path high_depth_regions_tbi

    output:
    tuple val(pdid), val(sample_id_ls), val(match_normal_id), path("tmpvaf_${sample_id_ls[0]}_${chrom}"), val(chrom)

    script:
    // meta = ['pdid':pdid, 'sample_id_ls':sample_id_ls, 'match_normal_id': match_normal_id]
    """
    rm -f "${pdid}_${mut_type}_vaf.tsv"
    rm -rf tmpvaf_*
    cgpVaf.pl -d . -o . -a ${mut_type} -g  ${fasta} -hdr  ${high_depth_regions} --vcf  ${vcf} --normal_bam ${bam_match} --tumour_bam ${bam} --normal_name ${match_normal_id} --tumour_name ${sample_id_ls.join(" ")} -chr ${chrom}
    mv "tmpvaf_${sample_id_ls[0]}" "tmpvaf_${sample_id_ls[0]}_${chrom}"
    """
}
