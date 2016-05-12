data {
  int<lower=1> N;  // Sample
  int<lower=0,upper=1> y[N];  // response
  vector[N] g;  // genotype
  real prior;
}

parameters {
  real p;  // variant effect
  real t;
  real<lower=0.01> sigma;
} 

transformed parameters {
  real mu;
  mu <- t * sigma;
}

model {
  p ~ normal(mu, sigma);
  t ~ normal(prior, 1);
  sigma ~ inv_gamma(2, 1);


  for (n in 1:N)
    y[n] ~ bernoulli_logit(p * g[n]); 
}

