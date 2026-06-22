# =============================================================================
# main_bdw_plugin.R  --  Bandwidth-comparison Monte Carlo (Table 3), plug-in
#
# Section 5 bandwidth study. Runs the Monte Carlo for the plug-in bandwidth
# selector (bw_plugin.R) across DGPs 1, 2, 3 (delta = c(0, 1, 0.5)[dgp]) at
# n = 1000, jobs 0-9, by sourcing funs/sim_bdw_plugin.R.
#
# Companion drivers main_bdw_loocv.R and main_bdw_rcv.R run the LOOCV and RCV
# selectors. Sources the shared core backend from ../core; results are written
# to results/plugin/ and aggregated into Table 3 by
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
seed1   <- 123       # Set initial seed (guaranteed reproducibility)
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
# Paths via the 'here' package (https://here.r-lib.org): resolved from the project
# root (the folder containing the '.here' sentinel / .git), so this script runs
# from ANY working directory with no editing. `address` is this scenario folder;
# CORE is the shared C++/R backend.
address <- here::here("simulation", "Bandwidth Choice (Table 3)")
CORE    <- here::here("simulation", "core")

#-----------------------------------------------------------------------------
# Source Rcpp file
Rcpp::sourceCpp(file.path(CORE, "np_lp.cpp"))

#-----------------------------------------------------------------------------
# Source R functions
source(file.path(CORE, "dgps_drdid.R"))
source(file.path(CORE, "locpol_funs.R"))
source(file.path(CORE, "dr_did.R"))
source(file.path(CORE, "bw_plugin.R"))

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Set seed
set.seed(seed1)

#-----------------------------------------------------------------------------
################################################################
# Run the simulations for all DGPs and all sample sizes
for (job in 0:9){
  for (nn in 3:3){
    for (dgp in 1:3){# 1: With Compositional Shifts; 2: No Shifts; 3: local alternative with delta = 0.5
        #set delta for local alternative with delta = 0.5
        delta <- c(0, 1, 0.5)[dgp]
        #set sample size
        if(nn==1) n <- 200
        if(nn==2) n <- 500
        if(nn==3) n <- 1000
        if(nn==4) n <- 5000
    
        # get system time
        ptm <- proc.time()
    
        #do the Monte Carlo
        source(file.path(address, "funs", "sim_bdw_plugin.R"))
    
        # return propensity score estimates
        time_used = proc.time() - ptm
    
        print(time_used)
    }
  }
}
