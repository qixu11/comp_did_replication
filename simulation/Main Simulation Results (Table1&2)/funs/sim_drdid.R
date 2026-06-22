# =============================================================================
# sim_drdid.R  --  Monte Carlo inner loop for the main simulation (Tables 1 & 2)
#
# Sourced by main_drdid.R with `dgp`, `n`, `job` (and `nrep`, `seed1`) in
# scope. Runs `nrep` replications: generates data via dgps_did(), fits the GPS
# and OR by cross-validated local polynomials, computes the stationary and
# non-stationary DR DiD estimators (ML/LS bandwidths), the TWFE benchmarks, and
# the Hausman-type test, storing all per-replication quantities in the
# 60-column matrix `mc`. Writes the per-job RData and summary CSVs to results/.
# =============================================================================
#########   MC Simulation file for Doubly Robust Diff-in-Diff       ########
############################################################################

tempData <- paste0(address, "/results", "/temp.dgp-",dgp,".n-",n,".job-", job,".RData")

pb <- txtProgressBar(max = nrep, style = 3)
progress <- function(n) setTxtProgressBar(pb,n)

#-----------------------------------------------------------------------------
#Start the MONTE CARLO loop
mc = matrix(0, nrep, 60)
for (imc in 1:nrep){
  #---------------------------------------------------------------------------
  #                          Data generating
  #---------------------------------------------------------------------------
  #Set seed
  iseed <- floor(seed1 + n*23 + job*1234 + imc)
  base::set.seed(iseed)
  
  #Generate the data
  data  = dgps_did(delta, n, xi_ps = 1)
  att.true = data$att.true
  att.unf  = data$att.unf
  att.unf.stnr  = data$att.unf.stnr
  ############################################################################
  ##############     DR ESTIMATOR with Stabilizing Weights          ##########
  ############################################################################
  list.control.ps = list(n_start = 1,
                         ps_min = 1e-5)
  
  
  ps.fit.ml = locpol.ps.fit(data, 
                            bws = NULL, 
                            cv_method= "cv.ml", 
                            lp_type= "logit",
                            list_control = list.control.ps)
  
  
  ps.fit.ls =  locpol.ps.fit(data, 
                             bws = NULL, 
                             cv_method= "cv.ls", 
                             lp_type= "logit",
                             list_control = list.control.ps)
  
  # Estimate OR functions
  list.control.or = list(n_start = 2,
                         bw_constrained = FALSE)
  
  or.fit  =  locpol.or.fit(data, 
                           bws = NULL, 
                           list_control = list.control.or)
  
  
  # Get the fitted values
  ps.hat.ml = ps.fit.ml$fitted.values 
  ps.hat.ls = ps.fit.ls$fitted.values 
  or.hat    = or.fit$fitted.values    
  
  
  ##############################################################################
  ######  Estimate ATTs
  ##############################################################################
  
  att_nst_ml <- att_st_ml <- att_nst_ls <- att_st_ls <- att_twfe1 <- att_twfe2  <- rep(0,2)
  se_nst_ml <- se_st_ml <- se_nst_ls <- se_st_ls <- se_twfe1 <- se_twfe2 <- rep(0,2)
  cp_nst_ml <- cp_st_ml <- cp_nst_ls <- cp_st_ls <- cp_twfe1 <- cp_twfe2 <- rep(0,2)
  len_nst_ml <- len_st_ml <- len_nst_ls <- len_st_ls <- len_twfe1 <- len_twfe2  <-  rep(0,2)
  inf.func.nst.ml <- inf.func.nst.ls <- inf.func.st.ml <- inf.func.st.ls <- matrix(0, n, 2)
  
  
  for (s in c(TRUE, FALSE)){
    
    i = 2 - as.numeric(s)
    
    ##############################################################################
    dr.nst.ml <- drdid_nonstnr(data$y, data$d, data$post, 
                               ps.hat.ml, or.hat,
                               stabilized = s, i.weights = NULL,
                               boot = FALSE, nboot = NULL,
                               inffunc = TRUE)
    
    # Get ATT and std. err. estimates
    att_nst_ml[i] <- dr.nst.ml$ATT
    se_nst_ml[i] <- dr.nst.ml$se
    # Whether the CI covers the true ATT (coverage probability)
    cp_nst_ml[i] <- as.numeric((dr.nst.ml$lci <= att.true) * (dr.nst.ml$uci >= att.true))
    # Length of confidence interval
    len_nst_ml[i] <- dr.nst.ml$uci - dr.nst.ml$lci
    # Infuence function
    inf.func.nst.ml[,i] <- dr.nst.ml$att.inf.func
    
    ##############################################################################
    dr.nst.ls <- drdid_nonstnr(data$y, data$d, data$post, 
                               ps.hat.ls, or.hat,
                               stabilized = s, i.weights = NULL,
                               boot = FALSE, nboot = NULL,
                               inffunc = TRUE)
    
    # Get ATT and std. err. estimates
    att_nst_ls[i] <- dr.nst.ls$ATT
    se_nst_ls[i] <- dr.nst.ls$se
    # Whether the CI covers the true ATT (coverage probability)
    cp_nst_ls[i] <- as.numeric((dr.nst.ls$lci <= att.true) * (dr.nst.ls$uci >= att.true))
    # Length of confidence interval
    len_nst_ls[i] <- dr.nst.ls$uci - dr.nst.ls$lci
    # Infuence function
    inf.func.nst.ls[,i] <- dr.nst.ls$att.inf.func
    
    ##############################################################################
    dr.st.ml <- drdid_stnr(data$y, data$d, data$post, 
                           ps.hat.ml, or.hat,
                           stabilized = s, i.weights = NULL,
                           boot = FALSE, nboot = NULL,
                           inffunc = TRUE)
    
    # Get ATT and std. err. estimates
    att_st_ml[i] <- dr.st.ml$ATT
    se_st_ml[i] <- dr.st.ml$se
    # Whether the CI covers the true ATT (coverage probability)
    cp_st_ml[i] <- as.numeric((dr.st.ml$lci <= att.true) * (dr.st.ml$uci >= att.true))
    #Length of confidence interval
    len_st_ml[i] <- dr.st.ml$uci - dr.st.ml$lci
    # Infuence function
    inf.func.st.ml[,i] <- dr.st.ml$att.inf.func
    
    ##############################################################################
    dr.st.ls <- drdid_stnr(data$y, data$d, data$post, 
                           ps.hat.ls, or.hat,
                           stabilized = s, i.weights = NULL,
                           boot = FALSE, nboot = NULL,
                           inffunc = TRUE)
    
    # Get ATT and std. err. estimates
    att_st_ls[i] <- dr.st.ls$ATT
    se_st_ls[i] <- dr.st.ls$se
    # Whether the CI covers the true ATT (coverage probability)
    cp_st_ls[i] <- as.numeric((dr.st.ls$lci <= att.true) * (dr.st.ls$uci >= att.true))
    #Length of confidence interval
    len_st_ls[i] <- dr.st.ls$uci - dr.st.ls$lci
    # Infuence function
    inf.func.st.ls[,i] <- dr.st.ls$att.inf.func
    
  }
  
  
  ##############################################################################
  ## Two-way fixed effect estimator -- baseline setup
  twfe1 <- twfe_did_rc(y = data$y, post = data$post, d= data$d, covariates = data$covariates)
  # Get ATT and std. err. estimates
  att_twfe1 <- twfe1$ATT
  se_twfe1  <- twfe1$se
  # Whether the CI covers the true ATT (coverage probability)
  cp_twfe1  <- as.numeric((twfe1$lci <= att.true) * (twfe1$uci >= att.true))
  # Length of confidence interval
  len_twfe1 <- twfe1$uci - twfe1$lci  
  
  ## Two-way fixed effect estimator -- saturated setup
  twfe2 <- twfe_did_rc(y = data$y, post = data$post, d= data$d, covariates = data$covariates_sat)
  # Get ATT and std. err. estimates
  att_twfe2 <- twfe2$ATT
  se_twfe2  <- twfe2$se
  # Whether the CI covers the true ATT (coverage probability)
  cp_twfe2  <- as.numeric((twfe2$lci <= att.true) * (twfe2$uci >= att.true))
  # Length of confidence interval
  len_twfe2 <- twfe2$uci - twfe2$lci  
  
  
  # Compute Wald statistics 
  test.ml.stb  <- n * (att_nst_ml[1] - att_st_ml[1])^2 / mean((inf.func.nst.ml[,1] - inf.func.st.ml[,1])^2)
  test.ls.stb  <- n * (att_nst_ls[1] - att_st_ls[1])^2 / mean((inf.func.nst.ls[,1] - inf.func.st.ls[,1])^2)
  test.ml.unstb  <- n * (att_nst_ml[2] - att_st_ml[2])^2 / mean((inf.func.nst.ml[,2] - inf.func.st.ml[,2])^2)
  test.ls.unstb  <- n * (att_nst_ls[2] - att_st_ls[2])^2 / mean((inf.func.nst.ls[,2] - inf.func.st.ls[,2])^2)
  
  
  
  # Get empirical rejection frequency 
  xi <- qchisq(c(0.9, 0.95, 0.99), df = 1)
  
  emp.rej.ml.stb <- test.ml.stb > xi
  emp.rej.ls.stb <- test.ls.stb > xi
  emp.rej.ml.unstb <- test.ml.unstb > xi
  emp.rej.ls.unstb <- test.ls.unstb > xi
  
  
  # Return output 
  mc[imc,] <- matrix(c(
    # true ATT
    att.true,
    # unfeasible ATT
    att.unf,
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
    test.ml.stb,    #36
    test.ls.stb,    #37
    test.ml.unstb,  #38
    test.ls.unstb,   #39
    
    # Empirical test size
    emp.rej.ml.stb, #40 - 42
    emp.rej.ls.stb, #43 - 45
    emp.rej.ml.unstb, #46 - 48
    emp.rej.ls.unstb, #49 - 51
    
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
    
    att.unf.stnr #60
  ), ncol = 60)
  
  progress(imc)
  
  # save temp results 
  if (imc%%100 == 0) {
    save.image(tempData)
  }
}

close(pb)



#-----------------------------------------------------------------------------
# Mean in the Monte Carlo
mean.mc <- base::colMeans(mc, na.rm = T)
# Median in MC
median.mc <- base::apply(mc, 2, FUN = median, na.rm=T)
# Bias
bias.mc <- base::colMeans(mc - mc[,1], na.rm = T)
median.bias.mc <- base::apply(mc - mc[,1], 2, FUN = median, na.rm=T)

# Standard deviation
sd.mc <- (base::colMeans(mc^2, na.rm = T) - base::colMeans(mc, na.rm = T)^2)^0.5

# RMSE
rmse.mc <- base::colMeans((mc - mc[1])^2, na.rm = T)^0.5

#MAE
mae.mc <- base::colMeans(abs(mc - mc[1]), na.rm = T)


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Create output tables
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#Summary Table with all results
# Create vector of summary statistics for all estimators
unf.nst.summary <- c(dgp,
                     n,
                     mean.mc[2], # average of the estimator
                     mean.mc[2] - mean.mc[1], # Bias
                     median.bias.mc[2],# Median Bias
                     rmse.mc[2], # RMSE
                     mae.mc[2], # MAE
                     sd.mc[2],  # Monte Carlo Std error
                     
                     NA, # Average Asy. variance
                     NA, # Empirical Coverage
                     NA, # Length of 95% Conf. Int.
                     NA  # Semiparametric efficiency Bound
                     
)

unf.st.summary <- c(dgp,
                    n,
                    mean.mc[52], #average of the estimator
                    mean.mc[52] - mean.mc[1], # Bias
                    median.bias.mc[52],#Median Bias
                    rmse.mc[52], #RMSE
                    mae.mc[52], #MAE
                    sd.mc[52],#Monte Carlo Std error
                    
                    NA, #Average Asy. variance
                    NA, # Empirical Coverage
                    NA, #Length of 95% Conf. Int.
                    NA # Semiparametric efficiency Bound
                    
)

# STABILIZED WEIGHTS
dr.nst.ml.summary.1 <-c(dgp,
                        n,
                        mean.mc[3], #average of the estimator
                        mean.mc[3] - mean.mc[1], # Bias
                        median.bias.mc[3],#Median Bias
                        rmse.mc[3], #RMSE
                        mae.mc[3], #MAE
                        sd.mc[3], #Monte Carlo Std error
                        
                        mean.mc[4],
                        mean.mc[5], # Empirical Coverage
                        mean.mc[6], #Length of 95% Conf. Int.
                        mean.mc[35] # Semiparametric efficiency Bound
)

dr.nst.ls.summary.1 <-    c(dgp,
                            n,
                            mean.mc[7], #average of the estimator
                            mean.mc[7] - mean.mc[1], # Bias
                            median.bias.mc[7],#Median Bias
                            rmse.mc[7], #RMSE
                            mae.mc[7], #MAE
                            sd.mc[7], #Monte Carlo Std error
                            
                            mean.mc[8],
                            mean.mc[9], # Empirical Coverage
                            mean.mc[10], #Length of 95% Conf. Int.
                            mean.mc[35] # Semiparametric efficiency Bound
)

dr.st.ml.summary.1 <-    c(dgp,
                           n,
                           mean.mc[11], #average of the estimator
                           mean.mc[11] - mean.mc[1], # Bias
                           median.bias.mc[11],#Median Bias
                           rmse.mc[11], #RMSE
                           mae.mc[11], #MAE
                           sd.mc[11], #Monte Carlo Std error
                           
                           mean.mc[12],
                           mean.mc[13], # Empirical Coverage
                           mean.mc[14], #Length of 95% Conf. Int.
                           mean.mc[35] # Semiparametric efficiency Bound
)

dr.st.ls.summary.1 <-    c(dgp,
                           n,
                           mean.mc[15], #average of the estimator
                           mean.mc[15] - mean.mc[1], # Bias
                           median.bias.mc[15],#Median Bias
                           rmse.mc[15], #RMSE
                           mae.mc[15], #MAE
                           sd.mc[15], #Monte Carlo Std error
                           
                           mean.mc[16],
                           mean.mc[17], # Empirical Coverage
                           mean.mc[18], #Length of 95% Conf. Int.
                           mean.mc[35] # Semiparametric efficiency Bound
)



# UNSTABILIZED WEIGHTS
dr.nst.ml.summary.2 <-    c(dgp,
                            n,
                            mean.mc[19], #average of the estimator
                            mean.mc[19] - mean.mc[1], # Bias
                            median.bias.mc[19],#Median Bias
                            rmse.mc[19], #RMSE
                            mae.mc[19], #MAE
                            sd.mc[19], #Monte Carlo Std error
                            
                            mean.mc[20],
                            mean.mc[21], # Empirical Coverage
                            mean.mc[22], #Length of 95% Conf. Int.
                            mean.mc[35] # Semiparametric efficiency Bound
)

dr.nst.ls.summary.2 <-    c(dgp,
                            n,
                            mean.mc[23], #average of the estimator
                            mean.mc[23] - mean.mc[1], # Bias
                            median.bias.mc[23],#Median Bias
                            rmse.mc[23], #RMSE
                            mae.mc[23], #MAE
                            sd.mc[23], #Monte Carlo Std error
                            
                            mean.mc[24],
                            mean.mc[25], # Empirical Coverage
                            mean.mc[26], #Length of 95% Conf. Int.
                            mean.mc[35] # Semiparametric efficiency Bound
)

dr.st.ml.summary.2 <-    c(dgp,
                           n,
                           mean.mc[27], #average of the estimator
                           mean.mc[27] - mean.mc[1], # Bias
                           median.bias.mc[27],#Median Bias
                           rmse.mc[27], #RMSE
                           mae.mc[27], #MAE
                           sd.mc[27], #Monte Carlo Std error
                           
                           mean.mc[28],
                           mean.mc[29], # Empirical Coverage
                           mean.mc[30], #Length of 95% Conf. Int.
                           mean.mc[35] # Semiparametric efficiency Bound
)

dr.st.ls.summary.2 <-    c(dgp,
                           n,
                           mean.mc[31], #average of the estimator
                           mean.mc[31] - mean.mc[1], # Bias
                           median.bias.mc[31],#Median Bias
                           rmse.mc[31], #RMSE
                           mae.mc[31], #MAE
                           sd.mc[31], #Monte Carlo Std error
                           
                           mean.mc[32],
                           mean.mc[33], # Empirical Coverage
                           mean.mc[34], #Length of 95% Conf. Int.
                           mean.mc[35] # Semiparametric efficiency Bound
)


# Two-way fixed effect
# baseline setting
twfe.base.summary      <- c(dgp,
                            n,
                            mean.mc[52], # average of the estimator
                            mean.mc[52] - mean.mc[1], # Bias
                            median.bias.mc[52],# Median Bias
                            rmse.mc[52], # RMSE
                            mae.mc[52],  # MAE
                            sd.mc[52],   # Monte Carlo Std error
                            
                            mean.mc[53], # Average Asy. variance
                            mean.mc[54], # Empirical Coverage
                            mean.mc[55], # Length of 95% Conf. Int.
                            mean.mc[35]  # Semiparametric efficiency Bound
)




# saturated setting
twfe.saturated.summary <-  c(dgp,
                             n,
                             mean.mc[56],  # average of the estimator
                             mean.mc[56] - mean.mc[1], # Bias
                             median.bias.mc[56], # Median Bias
                             rmse.mc[56],  # RMSE
                             mae.mc[56],   # MAE
                             sd.mc[56],    # Monte Carlo Std error
                             
                             mean.mc[57],  # Average Asy. variance
                             mean.mc[58],  # Empirical Coverage
                             mean.mc[59],  # Length of 95% Conf. Int.
                             mean.mc[35]   # Semiparametric efficiency Bound
)


# Summarize test results
# STABILIZED WEIGHTS
dr.test.ml.summary.1 <- c(dgp,
                          n,
                          "true",
                          mean.mc[36],
                          mean.mc[40],
                          mean.mc[41],
                          mean.mc[42])

dr.test.ls.summary.1 <- c(dgp,
                          n,
                          "true",
                          mean.mc[37],
                          mean.mc[43],
                          mean.mc[44],
                          mean.mc[45])

# UNSTABILIZED WEIGHTS
dr.test.ml.summary.2 <- c(dgp,
                          n,
                          "false",
                          mean.mc[38],
                          mean.mc[46],
                          mean.mc[47],
                          mean.mc[48])

dr.test.ls.summary.2 <- c(dgp,
                          n,
                          "false",
                          mean.mc[39],
                          mean.mc[49],
                          mean.mc[50],
                          mean.mc[51])


#-----------------------------------------------------------------------------
# Create the tables
scenarios =  c("unf.nstnr", "nstnr_ml", "nstnr_ls",
               "unf.stnr",  "stnr_ml",  "stnr_ls",
               "twfe_base", "twfe_saturated")
## STABILIZED WEIGHTS
mc.summary.1 <- rbind(unf.nst.summary,
                      dr.nst.ml.summary.1,
                      dr.nst.ls.summary.1,
                      unf.st.summary,
                      dr.st.ml.summary.1,
                      dr.st.ls.summary.1,
                      twfe.base.summary,
                      twfe.saturated.summary)


rownames(mc.summary.1) <- scenarios




colnames(mc.summary.1) <- c("DGP", "n", "Estimator", "Av.Bias", "Med.Bias","RMSE", "MAE", "MCSD",
                            "Asy.Var", "Coverage", "Lenth-CI", "Sem. Eff. Bound")


out1 <- paste0("mc.summary-dr-stb",".dgp-",dgp,".n-",n,".job-", job,".csv")
out1 <- paste0(address,"/results/", out1, sep="")

write.csv(mc.summary.1, file = out1)



## UNSTABILIZED WEIGHTS
mc.summary.2 <- rbind(unf.nst.summary,
                      dr.nst.ml.summary.2,
                      dr.nst.ls.summary.2,
                      unf.st.summary,
                      dr.st.ml.summary.2,
                      dr.st.ls.summary.2,
                      twfe.base.summary,
                      twfe.saturated.summary)


rownames(mc.summary.2) <- scenarios



colnames(mc.summary.2) <- c("DGP", "n", "Estimator", "Av.Bias", "Med.Bias","RMSE", "MAE", "MCSD",
                            "Asy.Var", "Coverage", "Lenth-CI", "Sem. Eff. Bound")


out2 <- paste0("mc.summary-dr-unstb",".dgp-",dgp,".n-",n,".job-", job,".csv")
out2 <- paste0(address,"/results/", out2, sep="")

write.csv(mc.summary.2, file = out2)


#  TEST RESULTS
mc.summary.3 <- rbind(dr.test.ml.summary.1,
                      dr.test.ls.summary.1,
                      dr.test.ml.summary.2,
                      dr.test.ls.summary.2)
colnames(mc.summary.3) <- c("DGP", "n", "Stabilized", "Av. Stats.", "Emp. Rej. Freq. (0.10)", "Emp. Rej. Freq. (0.05)", "Emp. Rej. Freq. (0.01)")
rownames(mc.summary.3) <- c("stb_ml", "stb_ls", "unstb_ml",  "unstb_ls")

out3 <- paste0("mc.summary-dr-test",".dgp-",dgp,".n-",n,".job-", job,".csv")
out3 <- paste0(address,"/results/", out3, sep="")

write.csv(mc.summary.3, file = out3)

#-----------------------------------------------------------------------------
# Save simulation results as RData

outData <- paste0(address, "/results", "/mc-dr.dgp-",dgp,".n-",n,".job-", job,".RData")
save.image(outData)
