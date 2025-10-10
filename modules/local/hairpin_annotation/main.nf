process hairpinAnnotation {
    label 'process_long'

    publishDir "${params.outdir}/filter_${mut_type}_out/${meta.pdid}", mode: params.publish_dir_mode

    // hairpin filter the mutation file
    input:
    tuple val(meta), path(bam), path(bai), path(bas), path(met), path(bam_match), path(bai_match), path(vcf), path(vcf_tbi)
    val mut_type
    path hairpin2_input_json
    path hairpin2_name_mapping

    output:
    tuple val(meta), path(vcf_annot_gz), emit: vcf_annot_gz
    tuple val(meta), path(vcf_annot_tbi), emit: vcf_annot_tbi

    script:
    vcf_annot = "${vcf.getName().tokenize(".").init().init().join(".")}.hairpin.vcf"
    vcf_annot_gz = "${vcf.getName().tokenize(".").init().init().join(".")}.hairpin.vcf.gz"
    vcf_annot_tbi = "${vcf.getName().tokenize(".").init().init().join(".")}.hairpin.vcf.gz.tbi"
    """
    hairpin2 --config ${hairpin2_input_json} --name-mapping ${hairpin2_name_mapping} ${vcf} ${bam} > ${vcf_annot}
    tabix_index.py --infile ${vcf_annot} --output_gz_path ${vcf_annot_gz} --preset vcf
    """
}
