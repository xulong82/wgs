//  Bayesian inference on ADSP project
//  ------------------------------------------------------------------------
//  Copyright: Xulong Wang <xulong.wang@jax.org>
//  License: GPLv3
//  ------------------------------------------------------------------------
//  STAN: hierarchical ordered logistic model

data {
  int<lower=1> N;  // Sample number
  int<lower=1> K;  // Ad categories
  int<lower=1> L;  // Family number
  int<lower=1,upper=K> Ad[N];  // Ad status
  int<lower=1,upper=L> Family[N];  // Family index
  int<lower=1,upper=2> Sex[N];  //  M:F
  vector[N] Age;  // Numerical
  vector[N] Snp;  // Numerical
}

transformed data { 
  vector[N] sexF;  // Relative to M

  for (i in 1:N) 
    sexF[i]  <- Sex[i]  == 2;
}

parameters {
  vector[L] alphaRaw;  // intercept
  real muAlpha;  // mu for intercept
  real<lower=0> sigmaAlpha;  // sigma for intercept

  real pAge;  // Age
  real pSex;  // Female
  real pSnp;  // SNP

  positive_ordered[K-2] cp1;  // cut points
} 

transformed parameters {
  vector[L] alpha;
  vector[K-1] cp2;  // cut points

  alpha <- muAlpha + sigmaAlpha * alphaRaw;

  cp2[1] <- 0;
  for (k in 2:(K-1)) 
    cp2[k] <- cp1[k-1];
}

model {
  vector[N] yhat;

  muAlpha ~ normal(0, 5);
  pAge ~ normal(0, 5);
  pSex ~ normal(0, 5);
  pSnp ~ normal(0, 5);

  alphaRaw ~ normal(0, 1);  // intercept

  for (n in 1:N) {
    yhat[n] <- alpha[Family[n]] + Age[n] * pAge + sexF[n] * pSex + Snp[n] * pSnp;
    Ad[n] ~ ordered_logistic(yhat[n], cp2);
  }
}

generated quantities {
  vector[N] yhat;
  vector[N] logLik;

  for (n in 1:N) {
    yhat[n] <- alpha[Family[n]] + Age[n] * pAge + sexF[n] * pSex + Snp[n] * pSnp;
    logLik[n] <- ordered_logistic_log(Ad[n], yhat[n], cp2);
  }
}

