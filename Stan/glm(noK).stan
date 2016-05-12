data {
  int<lower=1> N;  // Sample
  int<lower=1> K;  // Category
  int<lower=1> D;  // Predictor
  row_vector[D] cov[N];  // Covariate
  int<lower=1,upper=K> Ad[N];  // Response
}

parameters {
  ordered[K-1] c;
  vector[D] beta;
} 

model {
  c ~ normal(3, 1);
  beta ~ normal(0, 1);

  for (n in 1:N)
    Ad[n] ~ ordered_logistic(cov[n] * beta, c);
}

generated quantities {
  vector[N] lp;

  for (n in 1:N)
    lp[n] <- ordered_logistic_log(Ad[n], cov[n] * beta, c);
}

