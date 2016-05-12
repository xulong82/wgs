peak.detection <- function (x, y, cut) { # peak detection
  #--- parameters and variables
  x <- x * 1e-3  # Kb unit
  window = 1e3  # adjustable
  gap <- 1e3  # Inter peak distane: 5M 
  len <- length(x)
  fs <- rep(0, len)
  
  # --- Step 1: S function
  for (i in (window + 1):(len - window)) {
    maxL <- max(y[i] - y[i - window]:y[i-1])
    maxR <- max(y[i] - y[i + 1]:y[i + window])
    fs[i] <- (maxL + maxR) / 2
  }
  
  idx1 <- which(fs > 0)
  idx2 <- which(fs[idx1] - mean(fs[idx1]) > 1.5 * sd(fs[idx1]))
  
  # --- step 2: threshold
  idx3 <- which(y[idx1[idx2]] > cut)
  
  # --- step 3: local maximum
  index.cur <- idx1[idx2][idx3]
  index.cur2 <- index.cur1 <- index.cur
  Judge <- ifelse(length(index.cur) > 1, T, F)  # init
  
  while (Judge) {
    Judge <- F  # init
    for (i in 1:(length(index.cur1) - 1)) {
      if (x[index.cur1[i+1]] - x[index.cur1[i]] < gap) {
        Judge = T  # run again
        if (y[index.cur1[i+1]] >= y[index.cur1[i]]) {
          index.cur2[i] <- 0
        } else {
          index.cur2[i+1] <- 0
        }  # 2nd if
      }  # 1st if 
    }  # for loop
    index.cur <- index.cur2[index.cur2 != 0]
    index.cur2 <- index.cur1 <- index.cur
    Judge <- ifelse(Judge & (length(index.cur) > 1), T, F)
  }  # while control
  
  # --- head
  if (max(y[1:window]) > cut) index.cur <- c(which.max(y[1:window]), index.cur) 
  
  return(index.cur)
}

