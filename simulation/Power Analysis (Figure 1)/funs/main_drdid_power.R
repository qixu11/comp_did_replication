# =============================================================================
# main_drdid_power.R  --  Power-analysis Monte Carlo (Figure 1)
#
# Section 5 power study for the Hausman-type test of covariate stationarity.
# Fixes the DGP and sweeps the local-alternative mixing weight delta over
# seq(0, 1, 0.05) at n = 1000, running `nrep` replications per delta by
# sourcing funs/sim_power.R.
#
# Sources the shared core backend from ../core; per-delta results are written
# to results/ and turned into the power-curve plot (Figure 1) by
# simulation_summary_power.R.
# =============================================================================
#---------------------------------------------------
#      Stabalized-weighted ATT of DID
#      Repeated Cross Section
#      Doubly Robust DID estimators
#---------------------------------------------------
#-----------------------------------------------------------------------------
# Startup - clear memory, load packages, and set parameters
# Clear memory
rm(list=ls())
#-----------------------------------------------------------------------------
# Basic parameters for the simulation - Doesn't change over setups
seed1   <- 123       # Set initial seed (guaranteed reproducibility)
nrep    <- 500       # Monte Carlo replications
job     <- 0
#-----------------------------------------------------------------------------
# load the necessary libraries
library(foreach)
library(parallel)
library(doParallel)
library(doRNG)
library(Rcpp)
library(RcppArmadillo)
library(roptim)

#-----------------------------------------------------------------------------
# Set the Working Directory
address <- "./simulation/Power Analysis (Figure 1)"
setwd(address)

# Shared core backend (one copy for all simulation scenarios)
CORE <- file.path(dirname(address), "core")

#-----------------------------------------------------------------------------
# Source Rcpp file
Rcpp::sourceCpp(file.path(CORE, "np_lp.cpp"))

#-----------------------------------------------------------------------------
# Source R functions
source(file.path(CORE, "dgps_drdid.R"))
source(file.path(CORE, "locpol_funs.R"))
source(file.path(CORE, "dr_did.R"))

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Set seed
set.seed(seed1)

#-----------------------------------------------------------------------------
################################################################
# Run the simulations for all DGPs and all sample sizes
dgp = 3 # 1: With Compositional Shifts; 2: No Shifts; 3: local alternative with mixing weight delta^alt
for (nn in 3:3){
    #set delta for local alternative with mixing weight delta^alt
    for (delta in seq(0, 1, 0.05)){ #Controls deviation from the null of NCC condition
      #set sample size
      if(nn==1) n <- 200
      if(nn==2) n <- 500
      if(nn==3) n <- 1000
      if(nn==4) n <- 5000
  
      #get system time
      ptm <- proc.time()
  
      #do the Monte Carlo
      source("funs/sim_power.R")
  
      # return propensity score estimates
      time_used = proc.time() - ptm
  
      print(time_used)
  }
}
