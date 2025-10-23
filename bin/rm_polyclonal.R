#!/usr/bin/env Rscript

library(optparse)
library(logr)
options(mc.cores = parallel::detectCores())
options(stringsAsFactors = F)

# parsing arguments
option_list = list(
    make_option(c("--pdid"), action="store", default=NULL, type='character', help="sample ID"),
    make_option(c("--nr_path"), action="store", default=NULL, type='character', help="path to NR matrix (rows are variants, columns are samples)"),
    make_option(c("--nv_path"), action="store", default=NULL, type='character', help="path to NV matrix (rows are variants, columns are samples)"),
    make_option(c("--genotype_bin_path"), action="store", default=NULL, type='character', help="path to genotype_bin matrix (rows are variants, columns are samples)"),
    make_option(c("--clonality_path"), action="store", default=NULL, type='character', help="path to clonality results file that has one column for sample IDs and one column for clonality status"),
    make_option(c("--outdir"), action="store", default=NULL, type='character', help="path to output directory")
)
opt = parse_args(OptionParser(option_list=option_list, add_help_option=T))

pdid = opt$pdid
nr_path = opt$nr_path
nv_path = opt$nv_path
genotype_bin_path = opt$genotype_bin_path
clonality_path = opt$clonality_path
outdir = opt$outdir

log_open(file_name = paste0(outdir, '/', pdid, '.clonality.log'), logdir = F)
log_print(paste0('INPUTS: ', opt))

# read in data
nr = read.table(nr_path, header = T, check.names = F)
nv = read.table(nv_path, header = T, check.names = F)
genotype_bin = read.table(genotype_bin_path, header = T, check.names = F)
clonality = read.table(clonality_path, header = T, check.names = F)

clonal_samples = clonality$sample_id[clonality$clonality == 'clonal']

nr_clonal = nr[, colnames(nr) %in% clonal_samples, drop=F]
nv_clonal = nv[, colnames(nv) %in% clonal_samples, drop=F]
genotype_bin_clonal = genotype_bin[, colnames(genotype_bin) %in% clonal_samples, drop=F]

write.table(nr_clonal, file = paste0(outdir, '/NR_bbinom_filtered_clonal.txt'), sep = ' ', quote = F, row.names = T)
write.table(nv_clonal, file = paste0(outdir, '/NV_bbinom_filtered_clonal.txt'), sep = ' ', quote = F, row.names = T)
write.table(genotype_bin_clonal, file = paste0(outdir, '/genotype_bin_clonal.txt'), sep = ' ', quote = F, row.names = T)
