data {
  int<lower=1> N;  // Sample
  int<lower=1> D;  // Covariate
  row_vector[D] cov[N];  // Covariate
  vector[N] g;  // genotype
  int<lower=0,upper=1> Ad[N];  // phenotype
}

parameters {
  real a;  // Intercept
  real p;  // variant effect
  vector[D] beta;  // covariate effect
} 

model {
  a ~ normal(0, 1);
  p ~ normal(0, 1);
  beta ~ normal(0, 1);

  for (n in 1:N)
    Ad[n] ~ bernoulli_logit(a + p * g[n] + cov[n] * beta);
}

