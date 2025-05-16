**BSYNC-LCE**

Bayesian SYNchronization of climate records with annual Layer-Counting Errors

BSYNC-LCE is a bespoke version of BSYNC (Aquino-Lopez & Muschitiello, in prep.) developed for the automatic alignment of Greenland ice-core and speleothem climate records. It employs a probabilistic inversion model that incorporates prior knowledge of annual ice-layer counting uncertainty to estimate the age offset between the Greenland Ice Core Chronology 2005 (GICC05) and the speleothem U–Th timescale. Since the older portion of the IntCal20 radiocarbon calibration curve is largely based on the speleothem timescale, synchronizing Greenland ice-core and speleothem data with BSYNC-LCE facilitates the integration of ice-core, U–Th-dated, and radiocarbon-dated paleoclimate records onto a unified chronological framework spanning the past 50,000 years.

PLEASE CITE: 
1. Muschitiello, F. and Aquino-Lopez, M. A.: Continuous synchronization of the Greenland ice-core and U–Th timescales using probabilistic inversion, Clim. Past, 20, 1415–1435, https://doi.org/10.5194/cp-20-1415-2024, 2024.
2. Strawson, I., Faïn, X., Bauska, T.K., Muschitiello, F., Vladimirova, D.O, Tetzner, D.R., Humby, J., Thomas, E.R., Liu, P., Zhang, B., Grilli, R., Rhodes, R.H., Historical Southern Hemisphere biomass burning variability inferred from ice core carbon monoxide records, Proc. Natl. Acad. Sci. U.S.A. 121 (33) e2402868121, https://doi.org/10.1073/pnas.2402868121, 2024.

**What's in this Repository?**
- data : the folder contains 1. (input) ice core d18O data from GRIP on the GICC05 timescale; 2. (target) a stack of 14 U-Th-dated speleothem d18O records from Southeast Asia using a Monte Carlo Principal Component                   Method; 3. the GICC05 maximum counting error; 4. GICC05–Intcal match points from cosmogenic wiggle matching.
- BSYNC–LCE.R : is the main BSYNC script that will run the synchronization.
- load_data.R : load the data sets for synchronization (input, target, and other chronological information).
- functions.R : is a set of functions used by BSYNC–LCE.R. 
- PPL_functions.R : prior, proposal, and likelihood functions for the MCMC model used by BSYNC–LCE.R.
- gicc05_prep.R : prepare the GICC05 counting error priors used by BSYNC–LCE.R. 
- mcmc_setup.R : creates the Bayesian setup for the MCMC simulation using the library *BayesianTools*.
- posterior.R : extracts the MCMC results and estimates the posterior distributions of the model parameters.
- README.md : a brief description of BSYNC-LCE.
- LICENCE : is the BSYNC–LCE licence file.

