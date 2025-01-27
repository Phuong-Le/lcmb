process getChromCgpVaf {
    label 'process_tiny'

    input:
    path fai
    path high_depth_regions

    output:
    path chrom_list

    script:
    chrom_list = 'chrom_list.txt'
    """
    get_chrom_cgpvaf.py --fai ${fai} --hdr ${high_depth_regions} --outfile ${chrom_list}
    """
}
