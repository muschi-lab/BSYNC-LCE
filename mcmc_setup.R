##########################################################
# Set lower and upper prior boundaries 
lower <- c(nodes[1] - mce1, -rce)
upper <- c(nodes[1] + mce1, rce)

# Create the prior setup
prior <- createPrior(sampler = sampler, density = density, lower = lower, upper = upper, best = NULL)


# Collect prior, likelihood and posterior functions
bayesianSetup <- createBayesianSetup(likelihood = likelihood, prior = prior,
                                     catchDuplicates = FALSE)
# Setup MCMC simulation
settings <- list(iterations = bsync_settings$iters, burnin = bsync_settings$Burn, thin = bsync_settings$Thin) 
##########################################################
