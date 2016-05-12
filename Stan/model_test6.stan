data {
  int<lower=1> N;  // Sample
  int<lower=1> K;  // Category
  int<lower=1> D;  // Covariate
  row_vector[D] cov[N];  // Covariate
  cholesky_factor_cov[N] L;  // Covariance
  vector[N] g;  // genotype
  int<lower=1,upper=K> Ad[N];  // phenotype
  vector[2] prior;
}

parameters {
  real p;  // variant effect
  vector[D] beta;  // covariate effect
  ordered[K-1] c;  // cutpoints
  real t;
  real<lower=0.01> sigma;
  vector[N] z;
} 

transformed parameters {
  vector[N] u;
  u <- L * z;
}

model {
  c ~ normal(3, 1);
  beta ~ normal(0, 1);

  p ~ normal(mu, sigma);
  mu ~ normal(prior[1], 1);
  sigma ~ gamma(2, 1/prior[2]);

  z ~ normal(0, 1);

  for (n in 1:N)
    Ad[n] ~ ordered_logistic(p * g[n] + cov[n] * beta + u[n], c);
}

