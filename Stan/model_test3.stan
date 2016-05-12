data {
  int<lower=1> N;  // Sample
  int<lower=1> K;  // Category
  int<lower=1> D;  // Covariate
  row_vector[D] cov[N];  // Covariate
  cholesky_factor_cov[N] L;  // Covariance
  vector[N] g;  // genotype
  int<lower=1,upper=K> Ad[N];  // phenotype
  real prior;
}

parameters {
  real p;  // variant effect
  vector[D] beta;  // covariate effect
  ordered[K-1] c;  // cutpoints
  vector[N] z;
  real<lower=0.01> sigma;
} 

transformed parameters {
  vector[N] u;
  u <- L * z;
}

model {
  c ~ normal(3, 1);
  p ~ normal(sigma * prior, sigma);
  beta ~ normal(0, 1);
  z ~ normal(0, 1);
  sigma ~ normal(0, 1);

  for (n in 1:N)
    Ad[n] ~ ordered_logistic(p * g[n] + cov[n] * beta + u[n], c);
}

