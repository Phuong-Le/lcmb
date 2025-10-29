data {
  int<lower=1> K;          // number of mixture components
  int<lower=1> N;          // number of data points
  int<lower=0> t;          // truncated value
  array[N] int<lower=t> NR;         // read depth for each variant
  array[N] int NV;         // number of reads supporting variants for each variant
}

parameters {
  simplex[K] w;          // mixing proportions
  ordered[K] p;             // locations of mixture components
}

model {
  vector[K] log_w = log(w);  // cache log calculation
  p ~ uniform(0, 1);
  // p ~ beta(2, 5);
  for (n in 1:N) {
    vector[K] lps = log_w;
    for (k in 1:K) {
      lps[k] += binomial_lpmf(NV[n] | NR[n], p[k]);
      if (NR[n] < t) {
        lps[k] += negative_infinity();
      } else {
        // lps[k] += binomial_lpmf(NV[n] | NR[n], p[k]);
        lps[k] += -log_sum_exp(binomial_lpmf(t | NR[n], p[k]),
                              binomial_lccdf(t | NR[n], p[k]));
      }
    }
    target += log_sum_exp(lps);
  }
}

generated quantities {
  vector[K] log_w = log(w);  // cache log calculation
  vector[N] log_lik;
  for (n in 1:N) {
    vector[K] lps = log_w;
    for (k in 1:K) {
      lps[k] += binomial_lpmf(NV[n] | NR[n], p[k]);
      if (NR[n] < t) {
        lps[k] += negative_infinity();
      } else {
        // lps[k] += binomial_lpmf(NV[n] | NR[n], p[k]);
        lps[k] += -log_sum_exp(binomial_lpmf(t | NR[n], p[k]),
                              binomial_lccdf(t | NR[n], p[k]));
      }
    log_lik[n] = log_sum_exp(lps);
    }
  }
}




