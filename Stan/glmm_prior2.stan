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
  real t;  // prior
  real<lower=0.0001> sigma;
  vector[N] z;
} 

transformed parameters {
  real mu;
  vector[N] u;

  mu <- sigma * t;
  u <- L * z;
}

model {
  c ~ normal(3, 1);
  beta ~ normal(0, 1);

  sigma ~ inv_gamma(2, 1);

  t ~ normal(prior, 1);
  p ~ normal(mu, sigma);

  z ~ normal(0, 0.25);

  for (n in 1:N)
    Ad[n] ~ ordered_logistic(p * g[n] + cov[n] * beta + u[n], c);
}

