# WAIC  

Waic <- function(stanfit) {

  loglik <- extract (stanfit, "lp")$lp

  n1 <- nrow(loglik)  # number of samples
  n2 <- ncol(loglik)  # number of data points

  vars <- colSums((loglik - matrix(colMeans(loglik), n1, n2, byrow = T))^2) / (n1 - 1)
  pwaic <- sum(vars)  # effective parameter number
  lpd <- sum(log(colMeans(exp(loglik))))  # log pointwise predictive density
  waic <- -2 * (lpd - pwaic)

  return(list(WAIC = waic, Pwaic = pwaic, LPD = lpd))

}

