#################################################################################################################
# BSYNC-LCE: Bayesian SYNchronization with Layer-Counting Error priors
#################################################################################################################
# This is a bespoke version of BSYNC, Bayesian SYNchronization of proxy records (Aquino-Lopez & Muschitiello, in prep.), 
# for automatic alignment of Greenland ice-core and East Asian summer monsoon speleothem data. BSYNC-LCE leans on 
# a probabilistic inversion model that factors in prior knowledge of the annual ice-layer counting uncertainty and 
# estimates the age difference between the Greenland ice-core chronology 2005 (GICC05) and the U–Th timescale.
#################################################################################################################
# WRITTEN BY: F. Muschitiello (fm476@cam.ac.uk)

# PLEASE CITE: "Muschitiello, F. and Aquino-Lopez, M. A.: Continuous synchronization of the Greenland ice-core 
# and U–Th timescales using probabilistic inversion, Clim. Past, 20, 1415–1435, 
# https://doi.org/10.5194/cp-20-1415-2024, 2024." 
#################################################################################################################

###################################################
# Load libraries
require(BayesianTools)
require(truncnorm)
require(Rcpp)
require(rstan)
###################################################

###################################################
# Define working directory
user_dir <- "~/Desktop/BSYNC-LCE"
###################################################

###################################################
# Source script with user-defined functions
source(file = paste0(user_dir, "/", "functions.R"))
###################################################

###################################################
# Load: 
# 1. input (GRIP d18O) 
# 2. target (EASM pc1)
# 3. GICC05 maximum counting error
# 4. match points from cosmogenic wiggle matching
source(file = paste0(user_dir, "/", "load_data.R"))
###################################################

###################################################
# Settings for the inversion model and the MCMC
bsync_settings <- list(
  N = 200, # number of nodes (i.e. the input timescale is divided into N equally spaced segments)
  f = 5, # multiplying factor to allow 5x faster changes in RCE
  flex = 1.25, # multiplying factor to exploring outside MCE bounds (125%)
  inf = 2, # inflate by 2x the error of the target
  ta = 3, # t-distro parameter
  tb = 4, # t-distro parameter
  iters = 1.0e+6, # MCMC iterations
  Burn = 5000, # burn-in
  Thin = 5 # thinning
)
###################################################

###################################################
# Define number of segments along input age scale 
# and GICC05 counting error priors.
source(file = paste0(user_dir, "/", "gicc05_prep.R"))
###################################################

##################################################################################################
##################################################################################################

###################################################
# Proposal, Prior and Likelihood functions
source(file = paste0(user_dir, "/", "PPL_functions.R"))
###################################################

###################################################
# Setup MCMC
source(file = paste0(user_dir, "/", "mcmc_setup.R"))
###################################################

##########################################################
# RUN MCMC
chain <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DREAMzs",settings = settings)
##########################################################

##################################################################################################
##################################################################################################

##########################################################
# Extract MCMC samples
c <- getSample(chain, start = 0)
##########################################################

##########################################################
# Check convergence (for reliable results gD should be < 1.1)
print(Rhat(c))
##########################################################

###################################################
# Extract MCMC results and calculate posterior predictions
source(file = paste0(user_dir, "/", "posterior.R"))
###################################################

##################################################################################################
##################################################################################################

##########################################################
# Plot synchronized input
{
  par(mfrow=c(1,1))
  plot(align, type='l', xlim=c(10000,50000), 
       ylim=c(-1.25,1.25), xlab="U-Th years b2k", ylab="norm units")
  polygon(c(target[,1], rev(target[,1])), 
          c(target[,2]+(bsync_settings$inf * 2 * target[,3]), rev(target[,2] - (bsync_settings$inf * 2 * target[,3]))),
          col = adjustcolor("dodgerblue", alpha=0.05), border = NA)
  polygon(c(target[,1], rev(target[,1])), 
          c(target[,2] + (2 * target[,3]), rev(target[,2] - (2 * target[,3]))),
          col = adjustcolor("dodgerblue", alpha=0.15), border = NA)
  lines(target_med, col=4)
  abline(v=seq(10000,56000,2000), lty=3, lwd=0.5)
  legend("topright",c("GRIP (synchronized)", "EASM PC1"),
         lty=c(1, 1), col=c(1, 4), cex=0.75)
}
##########################################################

##########################################################
# Plot timescale offset ∆t  
##########################################################
{ 
  plot(quant[1,], quant[1,] - nodes, type='n', cex=0.25, pch=19, ylim=c(-750,750), xlim=c(10000,50000),
       xlab="U-Th years b2k", ylab="U-Th - GICC05 (years)")
  polygon(c(Dt$x, rev(Dt$x)), 
          c(Dt_95$y, rev(Dt_05$y)),
          col = adjustcolor("dodgerblue", alpha=0.15), border = NA)
  lines(Dt, lwd = 2, col = 4)
  lines(Dt_95, lwd = 2, lty = 3, col = 4)
  lines(Dt_05, lwd = 2, lty = 3, col = 4)
  
  # Cosmo match points
  points(TR$Age_yrs_b2k, TR$Dt, pch=19)
  arrows(TR$Age_yrs_b2k, TR$DT_95, 
         TR$Age_yrs_b2k, TR$DT_05,
         code=0)
  arrows(TR$Age_yrs_b2k - TR$age_error, TR$Dt, 
         TR$Age_yrs_b2k + TR$age_error, TR$Dt,
         code=0)

    # MCE bounds
  lines(GICC05_MCE[,1], GICC05_MCE[,2], lty=2, col = "grey")
  lines(GICC05_MCE[,1], -1*GICC05_MCE[,2], lty=2, col = "grey")
  
  abline(h=seq(-2000,2000,200), lty=3, lwd=0.5)
  abline(v=seq(10000,56000,2000), lty=3, lwd=0.5)
  
  legend("bottomright",c("∆t continuous", "cosmo match points (Muscheler et al. 2020)"),
         lty=c(1, NA), pch = c(NA, 19), col=c(4, 1), cex=0.75)
}
####################################################

##################################################################################################
##################################################################################################

####################################################
# Construct the output file path
output_file <- file.path(user_dir, "Dt.txt")
# Bind and round-up results
my_Dt <- cbind(round(Dt$x, 0), round(Dt$y, 0), round(Dt_95$y, 0), round(Dt_05$y, 0))

# Save to a CSV file
write.table(my_Dt, file = output_file, 
          row.names = FALSE, col.names = c("U-Th_years_b2k", "Dt_years", "crI_95", "crI_05"))
####################################################
