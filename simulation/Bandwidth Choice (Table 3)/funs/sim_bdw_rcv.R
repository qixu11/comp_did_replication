# =============================================================================
# sim_bdw_rcv.R  --  Monte Carlo inner loop for the bandwidth study (Table 3)
#
# Sourced by main_bdw_rcv.R with `dgp`, `n`, `job`, `delta` in scope. For each
# of `nrep` replications it generates data via dgps_did(), selects the
# local-polynomial bandwidths by RCV, fits the GPS/OR, computes the DR DiD
# estimators and the Hausman-type test, and stores them in the 60-column matrix
# `mc`. Writes per-job RData and summary CSVs to results/rcv/.
# The `_twfe` objects in this scenario are schema-preserving placeholders only:
# they keep the same output layout used by other simulation designs so the
# downstream summary code can read all scenarios uniformly.
# =============================================================================
# ---------------------------------------------------------------------------
# Monte Carlo simulation for doubly robust DID (RCV bandwidths)
# ---------------------------------------------------------------------------

temp_data <- paste0(address, "/results/rcv/", "temp_rcv.dgp-", dgp,
                   ".n-", n, ".job-", job, ".RData")

pb <- txtProgressBar(max = nrep, style = 3)
progress <- function(n) setTxtProgressBar(pb,n)

#-----------------------------------------------------------------------------
#Start the MONTE CARLO loop
mc <- matrix(0, nrep, 60)
for (imc in 1:nrep) {
  #---------------------------------------------------------------------------
  #                          Data generating
  #---------------------------------------------------------------------------
  #Set seed
  iseed <- floor(seed1 + n*23 + job*1234 + imc)
  base::set.seed(iseed)
  
  #Generate the data
  data <- dgps_did(delta, n, xi_ps = 1)
  att_true <- data$att.true
  att_unf <- data$att.unf
  att_unf_stnr <- data$att.unf.stnr
  ############################################################################
  ##############     DR ESTIMATOR with Stabilizing Weights          ##########
  ############################################################################
  list_control_ps <- list(n_start = 1,
                          ps_min = 1e-5,
                          cv_type = "rcv")

  
  ps_fit_ml <- locpol.ps.fit(data,
                             bws = NULL,
                             cv_method = "cv.ml",
                             lp_type = "logit",
                             list_control = list_control_ps)
  
  
  ps_fit_ls <- locpol.ps.fit(data,
                             bws = NULL,
                             cv_method = "cv.ls",
                             lp_type = "logit",
                             list_control = list_control_ps)
  
  # Estimate OR functions
  list_control_or <- list(n_start = 2,
                          bw_constrained = FALSE,
                          cv_type = "rcv")
  
  or_fit <- locpol.or.fit(data,
                          bws = NULL,
                          list_control = list_control_or)

  
  # Get the fitted values
  ps_hat_ml <- ps_fit_ml$fitted.values
  ps_hat_ls <- ps_fit_ls$fitted.values
  or_hat <- or_fit$fitted.values
  

  ##############################################################################
  ######  Estimate ATTs
  ##############################################################################
  
  att_nst_ml <- att_st_ml <- att_nst_ls <- att_st_ls <- att_twfe1 <- att_twfe2 <- rep(0, 2)
  se_nst_ml <- se_st_ml <- se_nst_ls <- se_st_ls <- se_twfe1 <- se_twfe2 <- rep(0, 2)
  cp_nst_ml <- cp_st_ml <- cp_nst_ls <- cp_st_ls <- cp_twfe1 <- cp_twfe2 <- rep(0, 2)
  len_nst_ml <- len_st_ml <- len_nst_ls <- len_st_ls <- len_twfe1 <- len_twfe2 <- rep(0, 2)
  inf_func_nst_ml <- inf_func_nst_ls <- inf_func_st_ml <- inf_func_st_ls <- matrix(0, n, 2)
  
  
  for (s in c(TRUE, FALSE)) {
    
    i = 2 - as.numeric(s)
    
    ##############################################################################
    dr_nst_ml <- drdid_nonstnr(data$y, data$d, data$post,
                               ps_hat_ml, or_hat,
                               stabilized = s, i.weights = NULL,
                               boot = FALSE, nboot = NULL,
                               inffunc = TRUE)
    
    # Get ATT and std. err. estimates
    att_nst_ml[i] <- dr_nst_ml$ATT
    se_nst_ml[i] <- dr_nst_ml$se
    # Whether the CI covers the true ATT (coverage probability)
    cp_nst_ml[i] <- as.numeric((dr_nst_ml$lci <= att_true) * (dr_nst_ml$uci >= att_true))
    # Length of confidence interval
    len_nst_ml[i] <- dr_nst_ml$uci - dr_nst_ml$lci
    # Infuence function
    inf_func_nst_ml[,i] <- dr_nst_ml$att.inf.func
    
    ##############################################################################
    dr_nst_ls <- drdid_nonstnr(data$y, data$d, data$post,
                               ps_hat_ls, or_hat,
                               stabilized = s, i.weights = NULL,
                               boot = FALSE, nboot = NULL,
                               inffunc = TRUE)
    
    # Get ATT and std. err. estimates
    att_nst_ls[i] <- dr_nst_ls$ATT
    se_nst_ls[i] <- dr_nst_ls$se
    # Whether the CI covers the true ATT (coverage probability)
    cp_nst_ls[i] <- as.numeric((dr_nst_ls$lci <= att_true) * (dr_nst_ls$uci >= att_true))
    # Length of confidence interval
    len_nst_ls[i] <- dr_nst_ls$uci - dr_nst_ls$lci
    # Infuence function
    inf_func_nst_ls[,i] <- dr_nst_ls$att.inf.func
    
    ##############################################################################
    dr_st_ml <- drdid_stnr(data$y, data$d, data$post,
                           ps_hat_ml, or_hat,
                           stabilized = s, i.weights = NULL,
                           boot = FALSE, nboot = NULL,
                           inffunc = TRUE)
    
    # Get ATT and std. err. estimates
    att_st_ml[i] <- dr_st_ml$ATT
    se_st_ml[i] <- dr_st_ml$se
    # Whether the CI covers the true ATT (coverage probability)
    cp_st_ml[i] <- as.numeric((dr_st_ml$lci <= att_true) * (dr_st_ml$uci >= att_true))
    #Length of confidence interval
    len_st_ml[i] <- dr_st_ml$uci - dr_st_ml$lci
    # Infuence function
    inf_func_st_ml[,i] <- dr_st_ml$att.inf.func
    
    ##############################################################################
    dr_st_ls <- drdid_stnr(data$y, data$d, data$post,
                           ps_hat_ls, or_hat,
                           stabilized = s, i.weights = NULL,
                           boot = FALSE, nboot = NULL,
                           inffunc = TRUE)
    
    # Get ATT and std. err. estimates
    att_st_ls[i] <- dr_st_ls$ATT
    se_st_ls[i] <- dr_st_ls$se
    # Whether the CI covers the true ATT (coverage probability)
    cp_st_ls[i] <- as.numeric((dr_st_ls$lci <= att_true) * (dr_st_ls$uci >= att_true))
    #Length of confidence interval
    len_st_ls[i] <- dr_st_ls$uci - dr_st_ls$lci
    # Infuence function
    inf_func_st_ls[,i] <- dr_st_ls$att.inf.func
    
  }
  
  
  ##############################################################################
  ## Two-way fixed effect estimator -- baseline setup  # Placeholders
  twfe1 <- 0 
  # Get ATT and std. err. estimates
  att_twfe1 <- 0
  se_twfe1  <- 0
  # Whether the CI covers the true ATT (coverage probability)
  cp_twfe1  <- 0 
  # Length of confidence interval
  len_twfe1 <- 0 
  
  ## Two-way fixed effect estimator -- saturated setup  # Placeholders
  twfe2 <- 0 
  # Get ATT and std. err. estimates
  att_twfe2 <- 0 
  se_twfe2  <- 0 
  # Whether the CI covers the true ATT (coverage probability)
  cp_twfe2  <- 0 
  # Length of confidence interval
  len_twfe2 <- 0
  
  
  # Compute Wald statistics 
  test_ml_stb <- n * (att_nst_ml[1] - att_st_ml[1])^2 /
    mean((inf_func_nst_ml[, 1] - inf_func_st_ml[, 1])^2)
  test_ls_stb <- n * (att_nst_ls[1] - att_st_ls[1])^2 /
    mean((inf_func_nst_ls[, 1] - inf_func_st_ls[, 1])^2)
  test_ml_unstb <- n * (att_nst_ml[2] - att_st_ml[2])^2 /
    mean((inf_func_nst_ml[, 2] - inf_func_st_ml[, 2])^2)
  test_ls_unstb <- n * (att_nst_ls[2] - att_st_ls[2])^2 /
    mean((inf_func_nst_ls[, 2] - inf_func_st_ls[, 2])^2)
  
  
  
  # Get empirical rejection frequency 
  xi <- qchisq(c(0.9, 0.95, 0.99), df = 1)
  
  emp_rej_ml_stb <- test_ml_stb > xi
  emp_rej_ls_stb <- test_ls_stb > xi
  emp_rej_ml_unstb <- test_ml_unstb > xi
  emp_rej_ls_unstb <- test_ls_unstb > xi
  
  
  # Return output 
  mc[imc, ] <- matrix(c(
    # true ATT
    att_true,
    # unfeasible ATT
    att_unf,
    ### CROSS-VALIDATED BANDWIDTH
    ## STABLIZED WEIGHTS
    # Nonstationary DR estimator - likelihood-based CV
    att_nst_ml[1],  #3
    (se_nst_ml[1] * sqrt(n))^2,   #4
    cp_nst_ml[1],      #5
    len_nst_ml[1],    #6
    
    # Nonstationary DR estimator - least-squares-based CV
    att_nst_ls[1],  #7
    (se_nst_ls[1] * sqrt(n))^2,   #8
    cp_nst_ls[1],      #9
    len_nst_ls[1],    #10
    
    # Stationary DR estimator - likelihood-based CV
    att_st_ml[1],  #11
    (se_st_ml[1] * sqrt(n))^2,   #12
    cp_st_ml[1],      #13
    len_st_ml[1],    #14
    
    # Stationary DR estimator - least-squares-based CV
    att_st_ls[1],  #15
    (se_st_ls[1] * sqrt(n))^2,   #16
    cp_st_ls[1],      #17
    len_st_ls[1],    #18
    
    ## UNSTABLIZED WEIGHTS
    # Nonstationary DR estimator - likelihood-based CV
    att_nst_ml[2],  #19
    (se_nst_ml[2] * sqrt(n))^2,   #20
    cp_nst_ml[2],      #21
    len_nst_ml[2],    #22
    
    # Nonstationary DR estimator - least-squares-based CV
    att_nst_ls[2],  #23
    (se_nst_ls[2] * sqrt(n))^2,   #24
    cp_nst_ls[2],      #25
    len_nst_ls[2],    #26
    
    # Stationary DR estimator - likelihood-based CV
    att_st_ml[2],  #27
    (se_st_ml[2] * sqrt(n))^2,   #28
    cp_st_ml[2],      #29
    len_st_ml[2],    #30
    
    # Stationary DR estimator - least-squares-based CV
    att_st_ls[2],  #31
    (se_st_ls[2] * sqrt(n))^2,   #32
    cp_st_ls[2],      #33
    len_st_ls[2],    #34
    
    
    # Semiparametric Efficiency bound
    data$eff, #35
    
    # Wald test stats. 
    test_ml_stb,    #36
    test_ls_stb,    #37
    test_ml_unstb,  #38
    test_ls_unstb,   #39
    
    # Empirical test size
    emp_rej_ml_stb, #40 - 42
    emp_rej_ls_stb, #43 - 45
    emp_rej_ml_unstb, #46 - 48
    emp_rej_ls_unstb, #49 - 51
    
    # Two-way fixed effect estimators 
    # baseline
    att_twfe1,  #52
    (se_twfe1 * sqrt(n))^2,  #53
    cp_twfe1,  #54
    len_twfe1,   #55
    # saturated 
    att_twfe2,  #56
    (se_twfe2 * sqrt(n))^2,  #57
    cp_twfe2,  #58
    len_twfe2,   #59
    
    att_unf_stnr #60
  ), ncol = 60)
  
  progress(imc)
  
  # save temp results 
  if (imc %% 100 == 0) {
    save.image(temp_data)
  }
}

close(pb)



#-----------------------------------------------------------------------------
# Mean in the Monte Carlo
mean_mc <- base::colMeans(mc, na.rm = T)
# Median in MC
median_mc <- base::apply(mc, 2, FUN = median, na.rm=T)
# Bias
bias_mc <- base::colMeans(mc - mc[,1], na.rm = T)
median_bias_mc <- base::apply(mc - mc[,1], 2, FUN = median, na.rm=T)

# Standard deviation
sd_mc <- (base::colMeans(mc^2, na.rm = T) - base::colMeans(mc, na.rm = T)^2)^0.5

# RMSE
rmse_mc <- base::colMeans((mc - mc[1])^2, na.rm = T)^0.5

#MAE
mae_mc <- base::colMeans(abs(mc - mc[1]), na.rm = T)


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Create output tables
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#Summary Table with all results
# Create vector of summary statistics for all estimators
unf_nst_summary <- c(dgp,
                     n,
                     mean_mc[2], # average of the estimator
                     mean_mc[2] - mean_mc[1], # Bias
                     median_bias_mc[2],# Median Bias
                     rmse_mc[2], # RMSE
                     mae_mc[2], # MAE
                     sd_mc[2],  # Monte Carlo Std error

                     NA, # Average Asy. variance
                     NA, # Empirical Coverage
                     NA, # Length of 95% Conf. Int.
                     NA  # Semiparametric efficiency Bound

)

unf_st_summary <- c(dgp,
                    n,
                    mean_mc[52], #average of the estimator
                    mean_mc[52] - mean_mc[1], # Bias
                    median_bias_mc[52],#Median Bias
                    rmse_mc[52], #RMSE
                    mae_mc[52], #MAE
                    sd_mc[52],#Monte Carlo Std error

                    NA, #Average Asy. variance
                    NA, # Empirical Coverage
                    NA, #Length of 95% Conf. Int.
                    NA # Semiparametric efficiency Bound

)

# STABILIZED WEIGHTS
dr_nst_ml_summary_1 <-c(dgp,
                        n,
                        mean_mc[3], #average of the estimator
                        mean_mc[3] - mean_mc[1], # Bias
                        median_bias_mc[3],#Median Bias
                        rmse_mc[3], #RMSE
                        mae_mc[3], #MAE
                        sd_mc[3], #Monte Carlo Std error

                        mean_mc[4],
                        mean_mc[5], # Empirical Coverage
                        mean_mc[6], #Length of 95% Conf. Int.
                        mean_mc[35] # Semiparametric efficiency Bound
)

dr_nst_ls_summary_1 <-    c(dgp,
                            n,
                            mean_mc[7], #average of the estimator
                            mean_mc[7] - mean_mc[1], # Bias
                            median_bias_mc[7],#Median Bias
                            rmse_mc[7], #RMSE
                            mae_mc[7], #MAE
                            sd_mc[7], #Monte Carlo Std error

                            mean_mc[8],
                            mean_mc[9], # Empirical Coverage
                            mean_mc[10], #Length of 95% Conf. Int.
                            mean_mc[35] # Semiparametric efficiency Bound
)

dr_st_ml_summary_1 <-    c(dgp,
                           n,
                           mean_mc[11], #average of the estimator
                           mean_mc[11] - mean_mc[1], # Bias
                           median_bias_mc[11],#Median Bias
                           rmse_mc[11], #RMSE
                           mae_mc[11], #MAE
                           sd_mc[11], #Monte Carlo Std error

                           mean_mc[12],
                           mean_mc[13], # Empirical Coverage
                           mean_mc[14], #Length of 95% Conf. Int.
                           mean_mc[35] # Semiparametric efficiency Bound
)

dr_st_ls_summary_1 <-    c(dgp,
                           n,
                           mean_mc[15], #average of the estimator
                           mean_mc[15] - mean_mc[1], # Bias
                           median_bias_mc[15],#Median Bias
                           rmse_mc[15], #RMSE
                           mae_mc[15], #MAE
                           sd_mc[15], #Monte Carlo Std error

                           mean_mc[16],
                           mean_mc[17], # Empirical Coverage
                           mean_mc[18], #Length of 95% Conf. Int.
                           mean_mc[35] # Semiparametric efficiency Bound
)



# UNSTABILIZED WEIGHTS
dr_nst_ml_summary_2 <-    c(dgp,
                            n,
                            mean_mc[19], #average of the estimator
                            mean_mc[19] - mean_mc[1], # Bias
                            median_bias_mc[19],#Median Bias
                            rmse_mc[19], #RMSE
                            mae_mc[19], #MAE
                            sd_mc[19], #Monte Carlo Std error

                            mean_mc[20],
                            mean_mc[21], # Empirical Coverage
                            mean_mc[22], #Length of 95% Conf. Int.
                            mean_mc[35] # Semiparametric efficiency Bound
)

dr_nst_ls_summary_2 <-    c(dgp,
                            n,
                            mean_mc[23], #average of the estimator
                            mean_mc[23] - mean_mc[1], # Bias
                            median_bias_mc[23],#Median Bias
                            rmse_mc[23], #RMSE
                            mae_mc[23], #MAE
                            sd_mc[23], #Monte Carlo Std error

                            mean_mc[24],
                            mean_mc[25], # Empirical Coverage
                            mean_mc[26], #Length of 95% Conf. Int.
                            mean_mc[35] # Semiparametric efficiency Bound
)

dr_st_ml_summary_2 <-    c(dgp,
                           n,
                           mean_mc[27], #average of the estimator
                           mean_mc[27] - mean_mc[1], # Bias
                           median_bias_mc[27],#Median Bias
                           rmse_mc[27], #RMSE
                           mae_mc[27], #MAE
                           sd_mc[27], #Monte Carlo Std error

                           mean_mc[28],
                           mean_mc[29], # Empirical Coverage
                           mean_mc[30], #Length of 95% Conf. Int.
                           mean_mc[35] # Semiparametric efficiency Bound
)

dr_st_ls_summary_2 <-    c(dgp,
                           n,
                           mean_mc[31], #average of the estimator
                           mean_mc[31] - mean_mc[1], # Bias
                           median_bias_mc[31],#Median Bias
                           rmse_mc[31], #RMSE
                           mae_mc[31], #MAE
                           sd_mc[31], #Monte Carlo Std error

                           mean_mc[32],
                           mean_mc[33], # Empirical Coverage
                           mean_mc[34], #Length of 95% Conf. Int.
                           mean_mc[35] # Semiparametric efficiency Bound
)


# Two-way fixed effect    # Placeholders
# baseline setting
twfe_base_summary      <- c(dgp,
                            n,
                            mean_mc[52], # average of the estimator
                            mean_mc[52] - mean_mc[1], # Bias
                            median_bias_mc[52],# Median Bias
                            rmse_mc[52], # RMSE
                            mae_mc[52],  # MAE
                            sd_mc[52],   # Monte Carlo Std error

                            mean_mc[53], # Average Asy. variance
                            mean_mc[54], # Empirical Coverage
                            mean_mc[55], # Length of 95% Conf. Int.
                            mean_mc[35]  # Semiparametric efficiency Bound
)




# saturated setting
twfe_saturated_summary <-  c(dgp,
                             n,
                             mean_mc[56],  # average of the estimator
                             mean_mc[56] - mean_mc[1], # Bias
                             median_bias_mc[56], # Median Bias
                             rmse_mc[56],  # RMSE
                             mae_mc[56],   # MAE
                             sd_mc[56],    # Monte Carlo Std error

                             mean_mc[57],  # Average Asy. variance
                             mean_mc[58],  # Empirical Coverage
                             mean_mc[59],  # Length of 95% Conf. Int.
                             mean_mc[35]   # Semiparametric efficiency Bound
)


# Summarize test results
# STABILIZED WEIGHTS
dr_test_ml_summary_1 <- c(dgp,
                          n,
                          "true",
                          mean_mc[36],
                          mean_mc[40],
                          mean_mc[41],
                          mean_mc[42])

dr_test_ls_summary_1 <- c(dgp,
                          n,
                          "true",
                          mean_mc[37],
                          mean_mc[43],
                          mean_mc[44],
                          mean_mc[45])

# UNSTABILIZED WEIGHTS
dr_test_ml_summary_2 <- c(dgp,
                          n,
                          "false",
                          mean_mc[38],
                          mean_mc[46],
                          mean_mc[47],
                          mean_mc[48])

dr_test_ls_summary_2 <- c(dgp,
                          n,
                          "false",
                          mean_mc[39],
                          mean_mc[49],
                          mean_mc[50],
                          mean_mc[51])


#-----------------------------------------------------------------------------
# Create the tables
scenarios =  c("unf.nstnr", "nstnr_ml", "nstnr_ls",
               "unf.stnr",  "stnr_ml",  "stnr_ls",
               "twfe_base", "twfe_saturated")
## STABILIZED WEIGHTS
mc_summary_1 <- rbind(unf_nst_summary,
                      dr_nst_ml_summary_1,
                      dr_nst_ls_summary_1,
                      unf_st_summary,
                      dr_st_ml_summary_1,
                      dr_st_ls_summary_1,
                      twfe_base_summary,     # Placeholders
                      twfe_saturated_summary # Placeholders
                      )


rownames(mc_summary_1) <- scenarios




colnames(mc_summary_1) <- c("DGP", "n", "Estimator", "Av.Bias", "Med.Bias","RMSE", "MAE", "MCSD",
                            "Asy.Var", "Coverage", "Lenth-CI", "Sem. Eff. Bound")


out1 <- paste0("mc.summary-dr-stb",".dgp-",dgp,".n-",n,".job-", job,".csv")
out1 <- paste0(address,"/results/rcv/", out1, sep="")

write.csv(mc_summary_1, file = out1)



## UNSTABILIZED WEIGHTS
mc_summary_2 <- rbind(unf_nst_summary,
                      dr_nst_ml_summary_2,
                      dr_nst_ls_summary_2,
                      unf_st_summary,
                      dr_st_ml_summary_2,
                      dr_st_ls_summary_2,
                      twfe_base_summary,     # Placeholders
                      twfe_saturated_summary # Placeholders
                      )


rownames(mc_summary_2) <- scenarios



colnames(mc_summary_2) <- c("DGP", "n", "Estimator", "Av.Bias", "Med.Bias","RMSE", "MAE", "MCSD",
                            "Asy.Var", "Coverage", "Lenth-CI", "Sem. Eff. Bound")


out2 <- paste0("mc.summary-dr-unstb",".dgp-",dgp,".n-",n,".job-", job,".csv")
out2 <- paste0(address,"/results/rcv/", out2, sep="")

write.csv(mc_summary_2, file = out2)


#  TEST RESULTS
mc_summary_3 <- rbind(dr_test_ml_summary_1,
                      dr_test_ls_summary_1,
                      dr_test_ml_summary_2,
                      dr_test_ls_summary_2)
colnames(mc_summary_3) <- c("DGP", "n", "Stabilized", "Av. Stats.", "Emp. Rej. Freq. (0.10)", "Emp. Rej. Freq. (0.05)", "Emp. Rej. Freq. (0.01)")
rownames(mc_summary_3) <- c("stb_ml", "stb_ls", "unstb_ml",  "unstb_ls")

out3 <- paste0("mc.summary-dr-test",".dgp-",dgp,".n-",n,".job-", job,".csv")
out3 <- paste0(address,"/results/rcv/", out3, sep="")

write.csv(mc_summary_3, file = out3)

#-----------------------------------------------------------------------------
# Save simulation results as RData

outData <- paste0(address, "/results/rcv/", "rlst_rcv.dgp-",dgp,".n-",n,".job-", job,".RData")
save.image(outData)
