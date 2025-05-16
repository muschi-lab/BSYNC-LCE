##########################################################
# Posterior age scale predictions
ages <- sapply(1:nrow(c), function(i){
  # extract and add up incremental dating steps
  d <- numeric(N + 1)
  d[1] <- c[i,1]
  for(j in 1:length(rce)){
    d[j+1] <- d[j] + interval + c[i,j+1]
  }
  result <- d
})
##########################################################

##########################################################
# Estimate age scale quantiles
quant <- apply(ages, 1, quantile, probs = c(0.5, 0.05, 0.95, 0.32, 0.68), na.rm =TRUE)
# Synchronized input
align <- cbind(approx(nodes, quant[1,], input[,1])$y, input[,2])
# Timescale offset ∆t 
Dt <- smooth.spline(quant[1,], quant[1,] - nodes, spar = 0.4)
Dt_95 <- smooth.spline(quant[1,], quant[2,] - nodes, spar = 0.4)
Dt_05 <- smooth.spline(quant[1,], quant[3,] - nodes, spar = 0.4)
##########################################################
