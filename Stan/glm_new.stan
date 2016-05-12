//  Copyright: Xulong Wang <xulong.wang@jax.org>
//  STAN: Generalized Linear Mixed Model

data {
  int<lower=1> N;  // Sample number
  int<lower=1> K;  // Ad categories
  matrix[N,N] KS;  // Kinship
  int<lower=1,upper=K> Ad[N];  // Ad status
  int<lower=1,upper=2> Sex[N];  //  M:F
  vector[N] Apoe2;  // Numerical
  vector[N] Apoe4;  // Numerical
  vector[N] Age;  // Numerical
  vector[N] Snp;  // Numerical
}

transformed data { 
  matrix[N,N] CK;  // Cholesky K

  vector[N] sexF;  // Female

  CK <- cholesky_decompose(KS);

  for (i in 1:N) {
    sexF[i]  <- Sex[i]  == 2;
  }
}

parameters {
  vector[N] z;  // ramdom effect
  real<lower=0.01, upper=1> sigmaU;

  simplex[K-1] cut1;  // cut points

  real alpha;  // Intercept
  real pAge;  // Age
  real pSex;  // Female
  real pSnp;  // SNP
  real pApoe2;  // APOE2
  real pApoe4;  // APOE4
} 

transformed parameters {
  vector[N] U;  // random effect
  ordered[K-1] cut2;
  real cutsum;

  cut2[1] <- 0;
  cutsum <- 0;

  for (cutp in 2:(K-1)) {
    cutsum <- cutsum + cut1[cutp-1];
    cut2[cutp] <- 10 * cutsum;
  }

  U <- sigmaU * (CK * z);
}

model {
  alpha ~ cauchy(0, 0.5);
  pAge ~ cauchy(0, 0.5);
  pSex ~ cauchy(0, 0.5);
  pSnp ~ cauchy(0, 0.5);
  pApoe2 ~ cauchy(0, 0.5);
  pApoe4 ~ cauchy(0, 0.5);
  sigmaU ~ cauchy(0, 0.5);

  cut1 ~ dirichlet(rep_vector(1, K-1));

  z ~ normal(0, 1);

  {
    vector[N] yhat;

    yhat <- U + alpha + Age*pAge + sexF*pSex + Snp*pSnp + Apoe2*pApoe2 + Apoe4*pApoe4;

    for (n in 1:N)
      Ad[n] ~ ordered_logistic(yhat[n], cut2);
  }
}

generated quantities {
  vector[N] logLik;

  {
    vector[N] yhat;

    yhat <- U + alpha + Age*pAge + sexF*pSex + Snp*pSnp + Apoe2*pApoe2 + Apoe4*pApoe4;
	
    for (n in 1:N) {
      logLik[n] <- ordered_logistic_log(Ad[n], yhat[n], cut2);
    }
  }
}

