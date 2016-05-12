//  Bayesian inference on ADSP project
//  Copyright: Xulong Wang <xulong.wang@jax.org>
//  ------------------------------------------------------------------------
//  STAN: hierarchical ordered logistic model

data {
  int<lower=1> N;  // Sample number
  int<lower=1> K;  // Ad categories
  matrix[N,N] KS;  // Kinship
  int<lower=1,upper=K> Ad[N];  // Ad status
  int<lower=1,upper=2> Sex[N];  //  M:F
  vector[N] Age;  // Numerical
  vector[N] Snp;  // Numerical
}

transformed data { 
  matrix[N,N] CK;  // Choleskty K
  vector[N] sexF;  // Relative to M

  CK <- cholesky_decompose(KS);

  for (i in 1:N)
    sexF[i]  <- Sex[i]  == 2;
}

parameters {
  vector[N] z;  // random effect
  real<lower=0> sigmaU;

  positive_ordered[K-1] cut;  // cut points

  real pAge;  // Age
  real pSex;  // Female
  real pSnp;  // SNP

} 

transformed parameters {
  vector[N] u;  // random effect

  for(i in 1:N) {
    real temp;
    temp <- 0;
    for(j in 1:i)
      temp <- temp + CK[i,j] * z[j];
    u[i] <- sigmaU * temp;
  }
}

model {
  pAge ~ cauchy(0, 1);
  pSex ~ cauchy(0, 1);
  pSnp ~ cauchy(0, 1);
  sigmaU ~ cauchy(0, 1);

  z ~ normal(0, 1);

  {
    vector[N] yhat;

    yhat <- u + Age * pAge + sexF * pSex + Snp * pSnp;

    for (n in 1:N)
      Ad[n] ~ ordered_logistic(yhat[n], cut);
  }
}

generated quantities {
  vector[N] logLik;

  {
    vector[N] yhat;
    yhat <- u + Age * pAge + sexF * pSex + Snp * pSnp;
 
    for (n in 1:N)
      logLik[n] <- ordered_logistic_log(Ad[n], yhat[n], cut);
  }
}
