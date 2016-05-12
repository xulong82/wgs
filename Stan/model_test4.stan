data {
  int<lower=1> N;  // Sample
  int<lower=1> K;  // Category
  int<lower=1> D;  // Covariate
  row_vector[D] cov[N];  // Covariate
  vector[N] g;  // genotype
  int<lower=1,upper=K> Ad[N];  // phenotype
  real prior;
}

parameters {
  vector[D] beta;  // covariate effect
  ordered[K-1] c;  // cutpoints

  real p;  // variant effect
  real<lower=0.001> scale;
} 

model {
  c ~ normal(3, 1);
  beta ~ normal(0, 1);

  scale ~ normal(0, 1);
  p ~ normal(prior * scale, scale);

  for (n in 1:N)
    Ad[n] ~ ordered_logistic(p * g[n] + cov[n] * beta, c);
}

