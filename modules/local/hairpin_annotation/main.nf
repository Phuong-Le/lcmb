process hairpinAnnotation {
    publishDir "${params.outdir}/filter_${mut_type}_out/${meta.pdid}", mode: params.publish_dir_mode

    // hairpin filter the mutation file
    input:
    tuple val(meta), path(bam), path(bai), path(bas), path(met), path(bam_match), path(bai_match), path(vcf), path(vcf_tbi)
    val hairpin_genome
    val mut_type

    output:
    tuple val(meta), path(vcf_annot)

    script:
    vcf_annot = "${vcf.getName().tokenize(".").init().init().join(".")}.hairpin.vcf"
    """
    hairpin -v ${vcf} -b ${bam} -o . -g ${hairpin_genome}
    """
}
