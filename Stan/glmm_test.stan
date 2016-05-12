data {
  int<lower=1> N;  // Sample
  int<lower=1> K;  // Category
  int<lower=1> D;  // Covariate
  row_vector[D] cov[N];  // Covariate
  cholesky_factor_cov[N] L;  // Covariance
  vector[N] g;  // genotype
  int<lower=1,upper=K> Ad[N];  // phenotype
}

parameters {
  real p;  // variant effect
  vector[D] beta;  // covariate effect
  ordered[K-1] c;  // cutpoints
  vector[N] z;
} 

model {
  vector[N] mu;
  for (n in 1:N) mu[n] <- 0;

  c ~ normal(3, 1);
  p ~ normal(0, 1);
  beta ~ normal(0, 1);
  z ~ multi_normal_cholesky(mu, L);

  for (n in 1:N)
    Ad[n] ~ ordered_logistic(p * g[n] + cov[n] * beta + z[n], c);
}

