process hairpinAnnotation {
    label 'process_long'

    publishDir "${params.outdir}/filter_${mut_type}_out/${meta.pdid}", mode: params.publish_dir_mode

    // hairpin filter the mutation file
    input:
    tuple val(meta), path(bam), path(bai), path(bas), path(met), path(bam_match), path(bai_match), path(vcf), path(vcf_tbi)
    val hairpin_genome
    val mut_type

    output:
    tuple val(meta), path(vcf_annot_gz)

    script:
    vcf_annot = "${vcf.getName().tokenize(".").init().init().join(".")}.hairpin.vcf"
    vcf_annot_gz = "${vcf.getName().tokenize(".").init().init().join(".")}.hairpin.vcf.gz"
    vcf_annot_tbi = "${vcf.getName().tokenize(".").init().init().join(".")}.hairpin.vcf.gz.tbi"
    """
    hairpin -v ${vcf} -b ${bam} -o . -g ${hairpin_genome}
    bgzip ${vcf_annot}
    tabix ${vcf_annot_gz}
    """
}
