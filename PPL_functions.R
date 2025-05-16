###################################################
# PROPOSAL
###################################################
# Uninformative proposal distribution
sampler <-  function(n=1){
  # First age point
  d1 <- runif(n, min = nodes[1] - mce1,
              max = nodes[1] + mce1)
  
  # RCE values
  drce <- numeric(N) 
  for(i in 1:length(rce)){
    drce[i] <- runif(n, min = -rce[i], 
                     max = rce[i])
  }
  
  d <- c(d1, drce)
  return(cbind(d))
}
##########################################################

###################################################
# PRIOR
###################################################
density <-  function(params){
  # PDF of the first age point
  d1 <- log(dtruncnorm(params[1], mean = nodes[1], sd = mce1/2,
                       a= nodes[1] - mce1, b = nodes[1] + mce1))
  
  # PDF of the rces
  drce <- numeric(N)
  for(i in 1:length(rce)){
    drce[i] <- dunif(params[i+1], min = -rce[i], max = rce[i], log = TRUE)
  }
  
  d <- c(d1, drce)
  return(sum(d))
}
##########################################################

##########################################################
# LIKELIHOOD
##########################################################
likelihood <- function(params){
  # Extract and add up incremental age steps
  d <- numeric(N + 1)
  d[1] <- params[1]
  for(i in 1:length(rce)){
    d[i+1] <- d[i] + interval + params[i+1]
  }
  
  # Apply simulated age scale to input record
  new_input <- approxima(nodes, d, input[,1])
  # Interpolate target and its stdev to resolution of simulated input
  new_target <- approxima_extrap(target[,1], target[,2], new_input)
  new_sd <- approxima_extrap(target[,1], target[,3], new_input)
  
  # Likelihood of observed values using t-distro
  likelihood <- tdistro(X = input[,2], Mu = new_target, sigma = bsync_settings$inf * new_sd, 
                        a = bsync_settings$ta, b = bsync_settings$tb)
  
  # Observed pre-15ky cosmo match points
  mp <- approxima(nodes, (d - nodes), TR[1:4,1])
  likelihood_cosmo <- dnorm(mp, mean = TR[1:4,3], sd = TR[1:4,6], log=TRUE) # conservative priors by setting the sdev to 2x the nominal value
  
  # sum the likelihood terms
  sumll = sum(likelihood, likelihood_cosmo)
  return(sumll)   
}
##########################################################