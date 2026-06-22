# =============================================================================
# main_bdw_loocv.R  --  Bandwidth-comparison Monte Carlo (Table 3), LOOCV
#
# Section 5 bandwidth study. Runs the Monte Carlo for the leave-one-out
# cross-validation (LOOCV) bandwidth selector across DGPs 1, 2, 3
# (delta = c(0, 1, 0.5)[dgp]) at n = 1000, jobs 0-9, by sourcing
# funs/sim_bdw_loocv.R.
#
# Companion drivers main_bdw_rcv.R and main_bdw_plugin.R run the RCV and
# plug-in selectors. Sources the shared core backend from ../core; results are
# written to results/loocv/ and aggregated into Table 3 by
# simulation_summary_table3.R.
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
nrep    <- 100       # Monte Carlo replications

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
# Working directory. Run this script from the repository ROOT: the relative default
# below is resolved to an absolute path automatically (no editing needed).
# Alternatively, set `address` to the absolute path of this scenario folder.
address <- normalizePath("./simulation/Bandwidth Choice (Table 3)")
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
################################################################
# Run the simulations for all DGPs and all sample sizes
for (nn in 3:3){
  for (dgp in 1:3){# 1: With Compositional Shifts; 2: No Shifts; 3: local alternative with delta = 0.5
    # set seed
    seed1   <- 123 * (dgp !=3) + 1234 * (dgp ==3)     # Set initial seed (guaranteed reproducibility)
    set.seed(seed1)
    
    for (job in 0:9){
      delta <- c(0, 1, 0.5)[dgp]
      #set sample size
      if(nn==1) n <- 200
      if(nn==2) n <- 500
      if(nn==3) n <- 1000
      if(nn==4) n <- 5000
      
      # get system time
      ptm <- proc.time()
      
      #do the Monte Carlo
      source("funs/sim_bdw_loocv.R")
      
      # return propensity score estimates
      time_used = proc.time() - ptm
      
      print(time_used)
    }
  }
}

