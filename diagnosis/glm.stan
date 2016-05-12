//  Bayesian inference on ADSP project
//  Copyright: Xulong Wang <xulong.wang@jax.org>
//  ------------------------------------------------------------------------
//  STAN: hierarchical ordered logistic model

data {
  int<lower=1> N;  // Sample number
  int<lower=1> K;  // Ad categories
  int<lower=1,upper=K> Ad[N];  // Ad status
  int<lower=1,upper=2> Sex[N];  //  M:F
  int<lower=1,upper=5> Apoe[N];  // 33:22:23:24:34
  vector[N] Age;  // Numerical
  vector[N] Snp;  // Numerical
  vector[N] U;  // Random effect w/o pU
}

transformed data { 
  vector[N] sexF;  // Relative to M
  vector[N] apoe22;  // Relative to 33
  vector[N] apoe23;
  vector[N] apoe24;
  vector[N] apoe34;

  for (i in 1:N) {
    sexF[i]  <- Sex[i]  == 2;
    apoe22[i] <- Apoe[i] == 2;
    apoe23[i] <- Apoe[i] == 3;
    apoe24[i] <- Apoe[i] == 4;
    apoe34[i] <- Apoe[i] == 5;
  }
}

parameters {
  positive_ordered[K-1] cut;  // cut points

  real pAge;  // Age
  real pSex;  // Female
  real pSnp;  // SNP
  real pApoe22;  // 22:33
  real pApoe23;  // 23:33
  real pApoe24;  // 24:33
  real pApoe34;  // 34:33
} 

model {
  pAge ~ cauchy(0, 1);
  pSex ~ cauchy(0, 1);
  pSnp ~ cauchy(0, 1);
  pApoe22 ~ cauchy(0, 1);
  pApoe23 ~ cauchy(0, 1);
  pApoe24 ~ cauchy(0, 1);
  pApoe34 ~ cauchy(0, 1);

  {
    vector[N] yhat;

    yhat <- U + Age * pAge + sexF * pSex + Snp * pSnp +
            apoe22 * pApoe22 + apoe23 * pApoe23 +
            apoe24 * pApoe24 + apoe34 * pApoe34;

    for (n in 1:N)
      Ad[n] ~ ordered_logistic(yhat[n], cut);
  }
}

generated quantities {
  vector[N] logLik;

  {
    vector[N] yhat;

    yhat <- U + Age * pAge + sexF * pSex + Snp * pSnp +
            apoe22 * pApoe22 + apoe23 * pApoe23 +
            apoe24 * pApoe24 + apoe34 * pApoe34;
 
    for (n in 1:N)
      logLik[n] <- ordered_logistic_log(Ad[n], yhat[n], cut);
  }
}

