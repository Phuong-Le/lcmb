#!/usr/bin/env bash

# parsing arguments
while [[ "$#" -gt 0 ]]
do
    case "$1" in
    --vcf )
        vcf="$2"
        shift 2
        ;;
    --vcfilter_config )
        vcfilter_config="$2"
        shift 2
        ;;
    --out_vcf_filtered )
        out_vcf_filtered="$2"
        shift 2
        ;;
    --)
        shift;
        break
        ;;
    *)
        echo "Unexpected option: $1"
        ;;
    esac
done

cmd="bcftools filter -i $(cat ${vcfilter_config}) -O z -o ${out_vcf_filtered} ${vcf}"
eval $cmd

# bgzip ${out_vcf_filtered}
# out_vcf_filtered_gz="${out_vcf_filtered}.gz"
# tabix ${out_vcf_filtered_gz}
bcftools index --tbi ${out_vcf_filtered}
