###################################################
## Load input and target records
# Records were: 
# 1. lowpass filtered with a cutoff period of 24 kyr to remove orbital-scale variability;
# 2. scaled between -1 and 1.

# GRIP d18O
input <- read.table(file = paste0(user_dir, "/data/", "GRIP_d18O_input.txt"), header = TRUE)
# South East Asian Monsoon speleo d18O PC1 
target <- read.table(file = paste0(user_dir, "/data/", "EASM_pc1_target.txt"), header = TRUE)
###################################################

###################################################
# GICC05 maximum counting error
GICC05_MCE <- read.table(file = paste0(user_dir, "/data/", "GICC05_mce.txt"), header = TRUE)
###################################################

###################################################
# Load match points from cosmogenic (14C-10Be) wiggle matching (Muscheler et al., 2020) 
TR <- read.table(file = paste0(user_dir, "/data/", "TR_match.txt"), header = TRUE)
###################################################
