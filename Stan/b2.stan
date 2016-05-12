data {
  int<lower=1> N;  // Sample
  int<lower=1> D;  // Covariate
  row_vector[D] cov[N];  // Covariate
  cholesky_factor_cov[N] L;  // Covariance
  vector[N] g;  // genotype
  int<lower=0,upper=1> Ad[N];  // phenotype
}

parameters {
  real a;  // Intercept
  real p;  // variant effect
  vector[D] beta;  // covariate effect
  real<lower=0.0001> sigma;
  vector[N] z;
} 

transformed parameters {
  vector[N] u;
  u <- L * z;
}


model {
  a ~ normal(0, 1);
  p ~ normal(0, 1);
  beta ~ normal(0, 1);
  sigma ~ inv_gamma(2, 1);
  z ~ normal(0, sigma);

  for (n in 1:N)
    Ad[n] ~ bernoulli_logit(a + p * g[n] + cov[n] * beta + u[n]);
}

