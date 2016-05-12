data {
  int<lower=1> N;  // Sample
  int<lower=1> D;  // Covariate
  row_vector[D] cov[N];  // Covariate
  vector[N] g;  // genotype
  vector[N] Ad;  // phenotype
}

parameters {
  real a;  // variant effect
  real p;  // variant effect
  vector[D] beta;  // covariate effect
  real<lower=0.01> sigma;
} 

model {
  a ~ normal(0, 1);
  beta ~ normal(0, 1);
  p ~ normal(0, 1);
  sigma ~ inv_gamma(2, 1);

  for (n in 1:N)
    Ad[n] ~ normal(a + p * g[n] + cov[n] * beta, sigma);
}

