data {
  int<lower=1> N;  // Sample
  int<lower=1> D;  // Covariate
  row_vector[D] cov[N];  // Covariate
  cholesky_factor_cov[N] L;  // Covariance
  vector[N] g;  // genotype
  int<lower=0,upper=1> Ad[N];  // phenotype
  vector[2] prior;
}

parameters {
  real p;  // variant effect
  vector[D] beta;  // covariate effect
  vector[N] z;
} 

transformed parameters {
  vector[N] u;
  u <- L * z;
}

model {
  p ~ normal(prior[1], prior[2]);
  beta ~ normal(0, 1);
  z ~ normal(0, 1);

  for (n in 1:N)
    Ad[n] ~ bernoulli_logit(p * g[n] + cov[n] * beta + u[n]);
}

generated quantities {
  vector[N] lp;

  for (n in 1:N)
    lp[n] <- bernoulli_logit_log(Ad[n], p * g[n] + cov[n] * beta + u[n]);
}

