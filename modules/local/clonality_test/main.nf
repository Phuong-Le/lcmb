process clonalityTest {
    label 'process_multicore'

    publishDir "${params.outdir}/filter_${mut_type}_out/${pdid}", pattern: "${sample_id}.clonality.log",  mode: params.publish_dir_mode
    publishDir "${params.outdir}/filter_${mut_type}_out/${pdid}", pattern: "${sample_id}_w.txt",  mode: params.publish_dir_mode
    publishDir "${params.outdir}/filter_${mut_type}_out/${pdid}", pattern: "${sample_id}_p.txt",  mode: params.publish_dir_mode
    publishDir "${params.outdir}/filter_${mut_type}_out/${pdid}", pattern: "${sample_id}_vaf_plot.pdf",  mode: params.publish_dir_mode


    input:
    tuple val(pdid), val(sample_id), path(vcf_path), path(nr_path), path(nv_path), path(genotype_bin_path)
    val min_good_reads
    val max_K
    val max_iter
    val nchains
    val clonal_threshold
    val proportion_pass_clonality
    val mut_type

    output:
    tuple val(pdid), path("${sample_id}_clonality.txt"), emit: clonality
    path "${sample_id}.clonality.log"
    path "${sample_id}_w.txt"
    path "${sample_id}_p.txt"
    path "${sample_id}_vaf_plot.pdf"


    script:
    """
    Rscript --vanilla ${projectDir}/bin/clonality_test.R --sample_id ${sample_id} --nr_path ${nr_path} --nv_path ${nv_path} --vcf_path ${vcf_path} --stan_path ${projectDir}/bin/truncated_binom_mixture.stan --truncated_value ${min_good_reads} --max_K ${max_K} --max_iter ${max_iter} --nchains ${nchains} --clonal_threshold ${clonal_threshold} --proportion_pass_clonality ${proportion_pass_clonality} --outdir .
    """
}
