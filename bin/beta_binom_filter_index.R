#!/usr/bin/env Rscript

options(stringsAsFactors = F)

# parsing arguments
{args = commandArgs(trailingOnly = TRUE)
# Function to parse named arguments
parse_arguments <- function(args) {
  parsed_args <- list()
  for (arg in args) {
    if (grepl("^--", arg)) {
      # Named argument found
      arg_split <- strsplit(arg, "=")[[1]]
      arg_name <- sub("^--", "", arg_split[1])
      arg_value <- arg_split[2]

      # Store the named argument in the list
      parsed_args[[arg_name]] <- arg_value
    }
  }
  return(parsed_args)
}
# Parse named arguments
parsed_args <- parse_arguments(args)
libpath = parsed_args$libpath
cgpvaf_out = parsed_args$cgpvaf_out
match_normal_id = parsed_args$match_normal_id
rho_threshold = parsed_args$rho_threshold
outdir = parsed_args$outdir
gender = parsed_args$gender
}


library(data.table)
library(VGAM)



source(paste0(libpath, '/Unmatched_NormSeq/germline_exact_binom.R'))
source(paste0(libpath, '/Unmatched_NormSeq/beta_binom_filter.R'))


vaf_data = read.table(cgpvaf_out, header=T)

ids = vaf_data[, c('Chrom', 'Pos', 'Ref', 'Alt')]

# If there is only one sample then no filtering is done
sample_vaf_cols = colnames(vaf_data)[grepl("VAF",colnames(vaf_data)) & colnames(vaf_data)!=paste0(match_normal_id,"_VAF")]
if (length(sample_vaf_cols) == 1) {

  no_filter = ids[, c('Chrom', 'Pos')]
  no_filter$Preceding_pos = no_filter$Pos - 1
  no_filter = no_filter[, c('Chrom', 'Preceding_pos', 'Pos')]
  write.table(no_filter, paste0(outdir,"/somatic_artefacts_one_sample_no_filter.bed"), quote=F, row.names=F , col.names=F , sep="\t")

} else if (length(sample_vaf_cols) < 1) {

  stop('there must be at least one sample')

} else {

  Muts = paste(vaf_data$Chrom,vaf_data$Pos,vaf_data$Ref,vaf_data$Alt,sep="_")
  Genotype = vaf_data[, grepl("VAF",colnames(vaf_data)) & colnames(vaf_data)!=paste0(match_normal_id,"_VAF")]
  NR = vaf_data[, grepl("DEP",colnames(vaf_data)) & colnames(vaf_data)!=paste0(match_normal_id,"_DEP")] # NR is the matrix with total depth (samples as columns, mutations rows
  NV = vaf_data[, grepl("MTR",colnames(vaf_data)) & colnames(vaf_data)!=paste0(match_normal_id,"_MTR")] # is matrix of reads supporting variants


  rownames(Genotype) = rownames(NV) = rownames(NR) = Muts

  samples = colnames(Genotype) = colnames(NR) = colnames(NV) = gsub("_VAF","",colnames(Genotype))

  XY_chromosomal = grepl("X|Y",Muts)
  autosomal = !XY_chromosomal
  xy_depth = mean(rowMeans(NR[XY_chromosomal,]))
  autosomal_depth = mean(rowMeans(NR[autosomal,]))

  if(is.null(gender)){
    gender='male'
    if(xy_depth>0.8*autosomal_depth) gender='female'
  }

  # Filter out germline using exact binomial test
  # germline=exact.binomial(gender=gender,NV=NV,NR=NR,cutoff = -5) #determine which variants are germline
  # # save indices to disk
  # germline_ids = ids[germline,]
  # somatic_ids = ids[!germline,]
  # write.table(germline_ids, paste0(outdir,"/germline_ids.txt"), row.names = F, col.names = T, quote=F)
  # write.table(somatic_ids, paste0(outdir,"/somatic_ids.txt"), row.names = F, col.names = T, quote=F)


  # Use beta-binomial filter on shared muts (unique muts would pass anyway) - filter out LCM artefacts
  NR_somatic=NR
  NV_somatic=NV
  somatic_ids = ids

  NR_somatic_nonzero=NR_somatic
  NR_somatic_nonzero[NR_somatic_nonzero==0]=1
  shared_muts=rowSums(NV_somatic>0)>1

  rho_est = beta.binom.filter(NR=NR_somatic_nonzero, NV=NV_somatic, shared_muts = shared_muts)
  flt_rho = rho_est < rho_threshold
  rho_df = data.frame(rho_est = rho_est, filter_out_by_rho = flt_rho)

  # save artefact filtered NR and NV to disk
  NR_somatic_noartefacts = NR_somatic_nonzero[!flt_rho,]
  NV_somatic_noartefacts = NV_somatic[!flt_rho,]
  write.table(NR_somatic_noartefacts, paste0(outdir,"/NR_bbinom_filtered.txt"), row.names = T, col.names = T, quote=F)
  write.table(NV_somatic_noartefacts, paste0(outdir,"/NV_bbinom_filtered.txt"), row.names = T, col.names = T, quote=F)

  #Convert genotype matrix in binary genotype
  XY_chromosomal=grepl("X|Y",rownames(NR_somatic_noartefacts))
  autosomal=!XY_chromosomal
  genotype_bin=as.matrix(NV_somatic_noartefacts/NR_somatic_noartefacts)
  if(gender=="male"){
    genotype_bin[autosomal,][genotype_bin[autosomal,]<0.1]=0
    genotype_bin[autosomal,][genotype_bin[autosomal,]>=0.3]=1
    genotype_bin[XY_chromosomal,][genotype_bin[XY_chromosomal,]<0.2]=0
    genotype_bin[XY_chromosomal,][genotype_bin[XY_chromosomal,]>=0.6]=1
    genotype_bin[genotype_bin>0&genotype_bin<1]=0.5
  }
  if(gender=="female"){
    genotype_bin[genotype_bin<0.1]=0
    genotype_bin[genotype_bin>=0.3]=1
    genotype_bin[genotype_bin>0&genotype_bin<1]=0.5
  }
  write.table(genotype_bin, paste0(outdir,"/genotype_bin.txt"), row.names = T, col.names = T, quote=F)



  # save rho and whether mutation is filtered to disk
  somatic_ids_rho = cbind(somatic_ids, rho_df)
  write.table(somatic_ids_rho, paste0(outdir,"/somatic_ids_rho.txt"), row.names = F, col.names = T, quote=F)

  # write ids for somatics and artefact filtered mutations
  somatic_artefacts_fltd = somatic_ids[!flt_rho, c('Chrom', 'Pos')]
  somatic_artefacts_fltd$Preceding_pos = somatic_artefacts_fltd$Pos - 1
  somatic_artefacts_fltd = somatic_artefacts_fltd[, c('Chrom', 'Preceding_pos', 'Pos')]

  write.table(somatic_artefacts_fltd, paste0(outdir,"/somatic_artefacts_filtered.bed"), quote=F, row.names=F , col.names=F , sep="\t")
}
