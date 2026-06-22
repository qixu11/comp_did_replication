###========================================================================
###========================================================================
### Doubly Robust DiD Estimates -- Sequeira (2016, AER)
###
### Replicates Panel B (DR DiD) and Panel C (Hausman test) of Table 4.
###
### Application: Sequeira (2016, AER) studies the effect of tariff
### reductions on bribe payments at the port of Maputo, Mozambique.
### Treatment: tariff_change_2008 (tariff reduced in 2008).
### Outcomes:  (1) Prob(bribe), (2) log(bribe), (3) log(bribe/shpt.val.),
###            (4) log(bribe/shpt.tonn.).
###
### The script:
### 1. Sources the local polynomial C++ backend and R wrapper functions.
### 2. Loops over four outcome variables.
### 3. For each outcome, estimates GPS and OR via LOOCV local polynomials.
### 4. Computes both the stationary (Sant'Anna & Zhao, 2020) and
###    non-stationary (Sant'Anna & Xu, 2026) DR DiD estimators.
### 5. Conducts the Hausman-type test for covariate stationarity.
### 6. Uses cluster bootstrap for inference (clustering on HS 4-digit codes).
### 7. Saves all results to results/result_sequeira_drdid.RData.
###========================================================================
###========================================================================

#-----------------------------------------------------------------------------
# Startup
#-----------------------------------------------------------------------------
rm(list = ls())

# Set the base directory to the application root.
# Users should modify this path to point to their local copy.
base_dir <- "./application" 

library("dplyr")
library("glmnet")

source(file.path(base_dir, "funs/pre_process.R"))
source(file.path(base_dir, "funs/locpol_funs.R"))
source(file.path(base_dir, "funs/dr_did.R"))

# Reproducibility
seed1 <- 1234
set.seed(seed1)

# Inference parameters
nboot <- 9999          # number of cluster bootstrap draws
alpha <- 0.95          # nominal coverage level
xi <- qchisq(alpha, df = 1)  # chi-squared critical value for Hausman test

# Import data (Sequeira, 2016 AER replication data)
data <- haven::read_dta(file.path(base_dir, "data/Bribes_Regression.dta"))

#-----------------------------------------------------------------------------
# Data preparation
#-----------------------------------------------------------------------------
# Fix a single data entry issue and recode variables
data$lba_value[2783] <- 0
data$monitor  <- data$monitor - 1   # recode to 0/1
data$psi      <- 2 - data$psi       # recode to 0/1
data$ltonnage <- log(data$tonnage + 1)

# Variable names
ynames <- c("bp", "lba", "lba_value", "lba_tonnage")
dname  <- "tariff_change_2008"  # treatment indicator
tname  <- "post_2008"           # post-period indicator
gname  <- "hc_4digits"          # cluster group (HS 4-digit code)

# Continuous covariates vary by outcome (different value/tonnage normalizations)
xcnames <- list(
  xcname1 = c("lvalue_tonnage", "tariff2007"),        # for Prob(bribe)
  xcname2 = c("lvalue_tonnage", "tariff2007"),        # for log(bribe)
  xcname3 = c("ltonnage", "tariff2007"),              # for log(bribe/shpt.val.)
  xcname4 = c("lvalue_shipment_metical", "tariff2007") # for log(bribe/shpt.tonn.)
)

# Unordered categorical covariates (common across all outcomes)
xunames <- c("differentiated", "agri", "perishable", "dfs",
             "day_w_arrival", "monitor", "psi", "rsa", "term")

# Select complete cases across all variables used in any specification
data <- data[stats::complete.cases(data[, c(ynames, tname, dname,
                                            "lvalue_tonnage", "ltonnage",
                                            "lvalue_shipment_metical",
                                            "tariff2007", xunames,
                                            gname, "clear_agent",
                                            "hc_group")]), ]

# Storage for results
ps.fit.ml <- or.fit <- sequeira_res <- list()

# Control parameters for bandwidth selection
list.control.ps <- list(n_start = 1,    # number of Nelder-Mead restarts
                        ps_min = 1e-2)  # propensity score trimming bound

list.control.or <- list(n_start = 1,         # number of restarts
                        bw_constrained = TRUE) # shared bandwidth across groups

###========================================================================
### Loop over four outcome variables
###========================================================================

for (iOutcome in 1:4) {
  # Construct the data list for outcome iOutcome
  dp <- pre_process_drdid(yname = ynames[iOutcome],
                          tname = tname,
                          dname = dname,
                          gname = gname,
                          xcnames = xcnames[[iOutcome]],
                          xunames = xunames,
                          data = data)
  

  # ---- Step 1: Estimate GPS (generalized propensity score) ----
  # Outcomes 1 and 2 share the same continuous covariates, so reuse GPS
  if (iOutcome != 2) {
    ps.fit.ml[[iOutcome]] <- locpol.ps.fit(dp,
                                           bws = NULL,
                                           cv_method = "cv.ml",
                                           lp_type = "logit",
                                           list_control = list.control.ps)
  } else {
    ps.fit.ml[[iOutcome]] <- ps.fit.ml[[1]]
  }
  fit.ps <- ps.fit.ml[[iOutcome]]$fitted.values

  # ---- Step 2: Estimate OR (outcome regression) ----
  or.fit[[iOutcome]] <- locpol.or.fit(dp,
                                      bws = NULL,
                                      list_control = list.control.or)
  fit.or <- or.fit[[iOutcome]]$fitted.values

  # ---- Step 3: DR DiD estimators ----

  # Stationary DR DiD (Sant'Anna & Zhao, 2020)
  dr_stnr <- drdid_stnr(
    dp$y,
    dp$d,
    dp$post,
    fit.ps,
    fit.or,
    stabilized = TRUE,
    boot = FALSE,
    inffunc = TRUE,
    alpha = alpha
  )
  # Non-stationary DR DiD (Sant'Anna & Xu, 2026)
  dr_nstnr <- drdid_nonstnr(
    dp$y,
    dp$d,
    dp$post,
    fit.ps,
    fit.or,
    stabilized = TRUE,
    boot = FALSE,
    inffunc = TRUE,
    alpha = alpha
  )
  
  # ---- Step 4: Hausman-type test for covariate stationarity ----
  # Wald statistic: n * (ATT_nstnr - ATT_stnr)^2 / Var(IF_nstnr - IF_stnr)
  wald_stat <- length(dp$d) * (dr_nstnr$ATT - dr_stnr$ATT)^2 /
               mean((dr_nstnr$att.inf.func - dr_stnr$att.inf.func)^2)
  stnr_test <- wald_stat > xi  # asymptotic rejection

  att_nstnr    <- dr_nstnr$ATT
  att_stnr     <- dr_stnr$ATT
  linrep_nstnr <- dr_nstnr$att.inf.func
  linrep_stnr  <- dr_stnr$att.inf.func

  # ---- Step 5: Cluster bootstrap inference ----
  # Jointly bootstrap both estimators for simultaneous inference
  linrep <- cbind(linrep_stnr, linrep_nstnr)
  dr.boot.list <- boot.drdid.cluster(linrep, dp$group, nboot)

  # Cluster-robust SE and CI for stationary estimator
  dr.boot  <- dr.boot.list$b.stnr
  sse_hat  <- IQR(dr.boot) / (qnorm(0.75) - qnorm(0.25))
  cv       <- quantile(abs(dr.boot / sse_hat), probs = alpha)
  se_stnr  <- sse_hat
  uci_stnr <- att_stnr + cv * sse_hat
  lci_stnr <- att_stnr - cv * sse_hat

  # Cluster-robust SE and CI for non-stationary estimator
  dr.boot   <- dr.boot.list$b.nstnr
  sse_hat   <- IQR(dr.boot) / (qnorm(0.75) - qnorm(0.25))
  cv        <- quantile(abs(dr.boot / sse_hat), probs = alpha)
  se_nstnr  <- sse_hat
  uci_nstnr <- att_nstnr + cv * sse_hat
  lci_nstnr <- att_nstnr - cv * sse_hat

  # Bootstrap test of covariate stationarity
  dr.boot.stat <- dr.boot.list$wald.stat
  p_value      <- mean(dr.boot.stat > wald_stat)
  stnr_btest   <- (p_value < 1 - alpha)
  
  # ---- Store results ----
  sequeira_res[[iOutcome]] <- list(
    att_nstnr  = att_nstnr,
    att_stnr   = att_stnr,
    se_nstnr   = dr_nstnr$se,
    lci_nstnr  = dr_nstnr$lci,
    uci_nstnr  = dr_nstnr$uci,
    se_stnr    = dr_stnr$se,
    lci_stnr   = dr_stnr$lci,
    uci_stnr   = dr_stnr$uci,
    wald_stat  = wald_stat,
    stnr_test  = stnr_test,
    sse_nstnr  = se_nstnr,
    slci_nstnr = lci_nstnr,
    suci_nstnr = uci_nstnr,
    sse_stnr   = se_stnr,
    slci_stnr  = lci_stnr,
    suci_stnr  = uci_stnr,
    p_value    = p_value,
    stnr_btest = stnr_btest
  )
  
}

# Save results
save.image(file.path(base_dir, "results/result_sequeira_drdid.RData"))