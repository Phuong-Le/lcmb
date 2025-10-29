#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggpubr)
library(logr)
library(rstan)
library(loo)
library(scales)
library(distributionsR)
library(optparse)
options(mc.cores = parallel::detectCores())
rstan_options(threads_per_chain = 1)
options(stringsAsFactors = F)

# parsing arguments
option_list = list(
  make_option(c("--sample_id"), action="store", default=NULL, type='character', help="sample ID"),
  make_option(c("--nr_path"), action="store", default=NULL, type='character', help="path to NR matrix (rows are variants, columns are samples)"),
  make_option(c("--nv_path"), action="store", default=NULL, type='character', help="path to NV matrix (rows are variants, columns are samples)"),
  make_option(c("--vcf_path"), action="store", default=NULL, type='character', help="Input VCF file"),
  make_option(c("--stan_path"), action="store", default=NULL, type='character', help="path to stan model file"),
  make_option(c("--truncated_value"), action="store", default=0, type='integer', help="truncated value of NV (NV >= truncated_value), default to 0"),
  make_option(c("--max_K"), action="store", default=3, type='integer', help="maximum number of clusters to trial"),
  make_option(c("--max_iter"), action="store", default=5, type='integer', help="maximum number of trials until modelling converges, default to 5"),
  make_option(c("--nchains"), action="store", default=5, type='integer', help="number of chains to run, default to 5"),
  make_option(c("--clonal_threshold"), action="store", default=0.25, type='numeric', help="threshold above which mutations are considered clonal"),
  make_option(c("--proportion_pass_clonality"), action="store", default=0.9, type='numeric', help="proportion of mutations required to pass clonality threshold"),
  make_option(c("--outdir"), action="store", default=NULL, type='character', help="path to output directory")
)
opt = parse_args(OptionParser(option_list=option_list, add_help_option=T))

sample_id = opt$sample_id
nr_path = opt$nr_path
nv_path = opt$nv_path
vcf_path = opt$vcf_path
stan_path = opt$stan_path
truncated_value = opt$truncated_value
max_K = opt$max_K # maximum
max_iter = opt$max_iter # maximum number of trials until modelling converges
nchains = opt$nchains
clonal_threshold = opt$clonal_threshold
proportion_pass_clonality = opt$proportion_pass_clonality
outdir = opt$outdir

log_open(file_name = paste0(outdir, '/', sample_id, '.clonality.log'), logdir = F)
log_print(paste0('INPUTS: ', opt))


# read in data
nr = read.table(nr_path, header = T, check.names = F)
nv = read.table(nv_path, header = T, check.names = F)
vcf = read.table(vcf_path, header = T, check.names = F)
mutations = paste(vcf[,1], vcf[,2], vcf[,4], vcf[,5], sep = '_') # 1, 2, 4, 5 correspond to chrom, pos, ref and pos
nr_sample = nr[mutations, sample_id]
nv_sample = nv[mutations, sample_id]


if (truncated_value < 0) stop(paste0(truncated_value, ' should be > 0'))
if (!all(nr_sample >= truncated_value)) log_warning(paste0(mutations[nr_sample < truncated_value], ' at ', nr_sample[nr_sample < truncated_value], ' depths and ', nv_sample[nr_sample < truncated_value], ' supporting reads have depths < ', truncated_value, '\n'))

if (!all(nv_sample >= truncated_value)) log_warning(paste0(mutations[nv_sample < truncated_value], ' at ', nr_sample[nv_sample < truncated_value], ' depths and ', nv_sample[nv_sample < truncated_value], ' supporting reads have number of supporting reads < ', truncated_value, '\n'))


log_print(paste0('minimum read depth: ', min(nr_sample)))

mutation_burden = length(nr_sample)
log_print(paste0('mutation burden: ', mutation_burden))
mean_depth = round(mean(nr_sample))
log_print(paste0('average read depth: ', mean_depth))

nv_sample_trunc_depth = nv_sample[which(nr_sample >= truncated_value)]
nr_sample_trunc_depth = nr_sample[which(nr_sample >= truncated_value)]
nv_sample_trunc_mtr = nv_sample_trunc_depth[which(nv_sample_trunc_depth >= truncated_value)]
nr_sample_trunc_mtr = nr_sample_trunc_depth[which(nv_sample_trunc_depth >= truncated_value)]


df = data.frame(NV = nv_sample_trunc_mtr,
                NR = nr_sample_trunc_mtr,
                vaf = nv_sample_trunc_mtr / nr_sample_trunc_mtr)
mutation_burden_trunc = nrow(df)

# apply mix model and model selection
# stan
mod_list = list()
for (K in 1:max_K) {
  log_print(paste0('trying ', K, ' cluster(s)'))
  stan_data = list(K = K,
                   N = mutation_burden_trunc,
                   t = truncated_value,
                   NR = df$NR,
                   NV = df$NV)
  # initial model fit
  fit_t <- stan(file = stan_path,
                 data = stan_data,
                 warmup = 2000,
                 iter = 3000,
                 chains = nchains,
                 cores = nchains)
  # has fit_t converged?
  mod_summary = summary(fit_t)
  rhat_max = max(mod_summary$summary[, 'Rhat'], na.rm = T)
  nparams = K*2 # 1 p and 1 w for each
  nchains = dim(mod_summary$c_summary)[3]
  log_print(paste0('number of initial successful chains: ', nchains))
  iter = 0
  while (rhat_max > 1.05 & iter < 5) {
    iter = iter + 1
    init_list = list()
    for (chain in 1:nchains) {
      init_list[[chain]] = as.list(mod_summary$c_summary[1:nparams, 1, chain]) # get the mean parameter estimate of each chain
    }
    fit_t <- stan(file = stan_path,
                   data = stan_data,
                   warmup = 1000,
                   iter = 2000,
                   chains = nchains,
                   cores = nchains,
                   init = init_list
    )
    mod_summary = summary(fit_t)
    nchains = dim(mod_summary$c_summary)[3]
    rhat_max = max(mod_summary$summary[, 'Rhat'], na.rm = T)
  }
  if (rhat_max > 1.05) {
    log_print('MODEL HAS NOT CONVERGED')
  } else {
    log_print(paste0('model convered after attempt ', iter, ', rhat_max = ', rhat_max))
  }
  mod_list[[K]] = fit_t
}

# model selection
loos_est = c()
loos_se = c()
for (K in 1:max_K) {
  mod = mod_list[[K]]
  summary(mod)
  posterior = extract(mod)
  if (is.null(posterior)) {
    log_print(paste0('failed modelling for '))
    loo_est = NA
    loo_se = NA
  } else {
    log_lik = extract_log_lik(mod)
    loo_mod = loo(log_lik)
    loo_est = loo_mod$estimates[1,1]
    loo_se = loo_mod$estimates[1,2]
  }
  loos_est = c(loos_est, loo_est)
  loos_se = c(loos_se, loo_se)
}

# rule of thumb: choosing the model with the lowest complexity whose estimated predictive performance is within one standard error of the best performance (https://dcl-model.stanford.edu/model_evaluation.html)
best_mod = which(loos_est == max(loos_est, na.rm = T))
log_print(paste0('Best model is K = ', best_mod))
loo_best = loos_est[best_mod]
loo_se_best = loos_se[which(loos_est == max(loos_est, na.rm = T))]
loo_threshold = loo_best + loo_se_best
chosen_K = which(loos_est < loo_threshold)[1]
chosen_mod = mod_list[[chosen_K]]
log_print(paste0('Chosen model is K = ', chosen_K, ' as it is within 1 standard error from the best model'))


posterior = extract(chosen_mod)
w_s = apply(posterior$w, MARGIN = 2, mean)
p_s = apply(posterior$p, MARGIN = 2, mean)

# so is the sample clonal?
clonal = ifelse(sum(w_s[which(p_s > clonal_threshold)]) > proportion_pass_clonality, 'clonal', 'polyclonal')
clonal_df = data.frame(sample_id = sample_id, clonality = clonal)
write.table(clonal_df, file = paste0(outdir, '/', sample_id, '_clonality.txt'), sep = '\t', row.names = F, quote = F, col.names = F)


# save model param estimates
w_df = as.data.frame(posterior$w)
colnames(w_df) = paste0('w', 1:length(w_s))
p_df = as.data.frame(posterior$p)
colnames(p_df) = paste0('p', 1:length(p_s))

write.table(w_df, file = paste0(outdir, '/', sample_id, '_w.txt'), sep = '\t', row.names = F, quote = F)
write.table(p_df, file = paste0(outdir, '/', sample_id, '_p.txt'), sep = '\t', row.names = F, quote = F)

# plotting
bw = 0.05
p = ggplot() +
  geom_histogram(data = df, mapping = aes(x = vaf), alpha = 0.7, position = 'identity', col = '#29335c', breaks = seq(0, 1, bw)) +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0,0))+
  ggtitle(paste0(sample_id, ', mean depth = ', mean_depth, ', snv burden = ', mutation_burden, ', truncated snv burden = ', mutation_burden_trunc)) +
  theme_pubr()

subtitles = paste0('p', 1:length(p_s), ' = ', round(p_s, 2), ', w', 1:length(w_s), ' = ', round(w_s, 2))
subtitles = paste(subtitles, collapse = '\n')

p = p + labs(subtitle = subtitles)

for (i in 1:length(p_s)) {
  w = w_s[i]
  peak = p_s[i]

  nmuts = round(mutation_burden_trunc*w)

  nr_sim = r_trunc_pois(n = nmuts, lambda = mean_depth, min_x = truncated_value)
  nv_sim = r_trunc_binom(n = nmuts, size = nr_sim, prob = peak, min_x = truncated_value)

  vaf = nv_sim / nr_sim
  vaf_val = hist(vaf, breaks = seq(0, 1, by = bw), xlim = c(0, 1))$mids
  vaf_freq = hist(vaf, breaks = seq(0, 1, by = bw), xlim = c(0,1))$counts
  vaf_df = data.frame(vaf = vaf_val, freq = vaf_freq)

  p = p +
    geom_smooth(data = vaf_df, aes(x = vaf, y = freq), se = F, col = hue_pal()(i)[i],  inherit.aes = F) + # use the default colours from ggplot
    scale_y_continuous(limits = c(0, NA), expand = c(0,0))
}

# save plot
ggsave(paste0(outdir, '/', sample_id, '_vaf_plot.pdf'), p)

log_close()


