#-------------------------------------------------
# Beta-binomial overdispersion - filtering out artefacts
# adapted from Tim Coorens - April 2018
# https://github.com/TimCoorens/Unmatched_NormSeq/tree/master
#-------------------------------------------------
require(VGAM)


estimateRho_gridml = function(NV_vec,NR_vec) {
  #Make sure depth is non-zero
  NV_vec=NV_vec[NR_vec>0]
  NR_vec=NR_vec[NR_vec>0]
  
  # Function to estimate maximum likelihood value of rho for beta-binomial
  rhovec = 10^seq(-6,-0.05,by=0.05) # rho will be bounded within 1e-6 and 0.89
  mu=sum(NV_vec)/sum(NR_vec)
  ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=NV_vec, size=NR_vec, rho=rhoj, prob=mu, log=T)))
  return(rhovec[ll==max(ll)][1])
}


beta.binom.filter = function(NR,NV,pval.cutoff=0.05,shared_muts=NULL){
  # Function to apply beta-binomial filter for artefacts. Works best on sets of
  # clonal samples (ideally >10 or so). As before, takes NV and NR as input. 
  # Optionally calculates pvalue of likelihood beta-binomial with estimated rho
  # fits better than binomial. This was supposed to protect against low-depth variants,
  # but use with caution. Returns logical vector with good variants = TRUE
  
  if (is.null(shared_muts)) shared_muts = rep(TRUE, nrow(NR))
  rho_est = pval = rep(NA,nrow(NR))
  for (k in 1:nrow(NR)){
    if (shared_muts[k] == FALSE) 
      {rho_est[k] = 0.89} else
      {rho_est[k]=estimateRho_gridml(NV_vec = as.numeric(NV[k,]),
                                      NR_vec = as.numeric(NR[k,]))}
    if (k%%1000==0){
      print(paste0(k, 'th mutation'))
    }
  }
  return(rho_est)
}