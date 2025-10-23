process concatClonality {
    label 'process_tiny'

    publishDir "${params.outdir}/filter_${mut_type}_out/${pdid}",  mode: params.publish_dir_mode


    input:
    tuple val(pdid), path(nr_path), path(nv_path), path(genotype_bin_path), path(clonality)
    val mut_type

    output:
    tuple val(pdid), path("${pdid}_clonality.txt"), emit: clonality_input


    script:
    """
    printf "sample_id\\tclonality\\n" > ${pdid}_clonality.txt
    cat ${clonality} >> ${pdid}_clonality.txt
    """
}
