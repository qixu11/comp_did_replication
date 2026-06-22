# =============================================================================
# main_drdid.R  --  Driver for the main Monte Carlo study (Tables 1 & 2)
#
# Section 5 main simulation. For each DGP (1 = non-stationary / compositional
# change, 2 = stationary) at n = 1000, runs the Monte Carlo over jobs 0-9
# (each job = `nrep` replications) by sourcing funs/sim_drdid.R.
#
# Every replication estimates the nonparametric DR DiD estimators (stationary
# and non-stationary, each with likelihood- and least-squares-cross-validated
# local-polynomial bandwidths), the two TWFE benchmarks, and the Hausman-type
# test for covariate stationarity.
#
# Sources the shared core backend (C++ local-polynomial routines + DGP and
# estimator wrappers) from ../core. Per-job results/summaries are written to
# results/; simulation_summary_table1&2.R aggregates them into Tables 1 & 2.
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
seed1   <- 123      # Set initial seed (guaranteed reproducibility)
nrep    <- 500       # Monte Carlo replications

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
address <- "./simulation/Main Simulation Results (Table1&2)"  # EDIT to your local simulation/<scenario> path
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
for (nn in 3:3){
  for (dgp in 1:2){# 1: With Compositional Shifts; 2: No Shifts
    #set delta for local alternative with delta = 0.5
    for (job in 0:9){
      delta <- c(0, 1)[dgp]
      #set sample size
      if(nn==1) n <- 200
      if(nn==2) n <- 500
      if(nn==3) n <- 1000
      if(nn==4) n <- 5000
      
      # get system time
      ptm <- proc.time()
      
      #do the Monte Carlo
      source("funs/sim_drdid.R")
      
      # return propensity score estimates
      time_used = proc.time() - ptm
      
      print(time_used)
    }
  }
}

