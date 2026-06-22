#-----------------------------------------------------------------------------
# DiD estimators and bootstrap inference for the application.
#
# This file provides:
#   drdid_nonstnr()       - Non-stationary DR DiD (Sant'Anna & Xu, 2023)
#   drdid_stnr()          - Stationary DR DiD (Sant'Anna & Zhao, 2020)
#   twfe_did_rc()         - Two-way fixed effects DiD (repeated cross-sections)
#   mboot.did()           - Mammen (1993) multiplier bootstrap
#   boot.did.cluster()    - Cluster-level multiplier bootstrap (single estimator)
#   boot.drdid.cluster()  - Cluster-level multiplier bootstrap (joint: stnr + nstnr)
#
# Note: these are application-specific versions. The simulation code in
# /funs/dr_did.R uses a slightly different interface (drdid_nonstationary).
#-----------------------------------------------------------------------------

#' Non-stationary DR DiD estimator (Sant'Anna & Xu, 2023)
#'
#' Uses four-group GPS weights to allow compositional changes across periods.
#'
#' @param y       Outcome vector.
#' @param d       Treatment indicator (1 = treated).
#' @param post    Post-period indicator (1 = post).
#' @param fit.ps  n x 4 GPS matrix (treat-post, treat-pre, ctrl-post, ctrl-pre).
#' @param fit.or  n x 4 OR fitted values (same column order).
#' @param stabilized Logical; use stabilized weights.
#' @param i.weights  Optional observation weights.
#' @param boot    Logical; use multiplier bootstrap.
#' @param nboot   Number of bootstrap draws.
#' @param inffunc Logical; return influence function.
#' @param alpha   Nominal coverage level.
#' @return List with ATT, se, uci, lci, boots, att.inf.func.
drdid_nonstnr <- function(y, d, post, fit.ps, fit.or,
                          stabilized = TRUE, i.weights = NULL,
                          boot = FALSE, nboot = NULL,
                          inffunc = FALSE, alpha = 0.95) {
  # D as vector
  d <- as.vector(d)
  # Sample size
  n <- length(d)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  
  
  #summary(fit.ps)
  pmin <- 0.0001
  fit.ps[fit.ps < pmin] = pmin
  fit.ps[fit.ps > (1-pmin)] = 1-pmin
  
  
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  
  # Get the PS estimates
  ps.fit.treat.post = fit.ps[,1]
  ps.fit.treat.pre  = fit.ps[,2]
  ps.fit.cont.post  = fit.ps[,3]
  ps.fit.cont.pre   = fit.ps[,4]
  #ps.fit.treat = ps.fit.treat.post + ps.fit.treat.pre
  
  # Get the OR estimates
  out.y.treat.post = fit.or[,1]
  out.y.treat.pre  = fit.or[,2]
  out.y.cont.post  = fit.or[,3]
  out.y.cont.pre   = fit.or[,4]
  
  #-----------------------------------------------------------------------------
  # Compute the ATT
  
  # First, the weights
  w.treat.post <-  d * post
  w.treat.pre  <-  d * (1 - post) * ps.fit.treat.post / ps.fit.treat.pre
  w.cont.post  <-  (1 - d) * post * ps.fit.treat.post / ps.fit.cont.post
  w.cont.pre   <-  (1 - d) * (1 - post) * ps.fit.treat.post / ps.fit.cont.pre
  
  
  
  # Elements of the influence function (summands)
  if (stabilized == TRUE)
  {
    eta.treat.post <- w.treat.post * (y - out.y.treat.post) / mean(w.treat.post)
    eta.treat.pre  <- w.treat.pre * (y - out.y.treat.pre) / mean(w.treat.pre)
    eta.cont.post  <- w.cont.post * (y - out.y.cont.post) / mean(w.cont.post)
    eta.cont.pre   <- w.cont.pre  * (y - out.y.cont.pre) / mean(w.cont.pre)
    eta.tau.x  <- w.treat.post * (out.y.treat.post -  out.y.treat.pre -
                                    out.y.cont.post  +  out.y.cont.pre) / mean(w.treat.post)
  } else 
  {
    eta.treat.post <- w.treat.post * (y - out.y.treat.post) / mean(w.treat.post)
    eta.treat.pre  <- w.treat.pre * (y - out.y.treat.pre) / mean(w.treat.post)
    eta.cont.post  <- w.cont.post * (y - out.y.cont.post) / mean(w.treat.post)
    eta.cont.pre   <- w.cont.pre  * (y - out.y.cont.pre) / mean(w.treat.post)
    eta.tau.x  <- w.treat.post * (out.y.treat.post -  out.y.treat.pre -
                                    out.y.cont.post  +  out.y.cont.pre) / mean(w.treat.post)
  }
  
  
  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)
  att.tau.x <-  mean(eta.tau.x)
  # estimate ATT 
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre)  + att.tau.x
  
  #-----------------------------------------------------------------------------
  # get the influence function and compute standard error
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # Influence function for the treated component
  inf.treat <- eta.treat.post - eta.treat.pre
  #-----------------------------------------------------------------------------
  # Influence function for the control component
  inf.cont <- eta.cont.post - eta.cont.pre
  #-----------------------------------------------------------------------------
  #get the influence function of the DR estimator (put all pieces together)
  dr.att.inf.func1 <- inf.treat - inf.cont
  #-----------------------------------------------------------------------------
  # Now, we only need to get the influence function of the adjustment terms
  # First, the terms as if all OR parameters were known
  inf.eff <- eta.tau.x - w.treat.post * dr.att / mean(w.treat.post)
  #-----------------------------------------------------------------------------
  #get the influence function of the locally efficient DR estimator (put all pieces together)
  dr.att.inf.func <- dr.att.inf.func1 + inf.eff
  
  #-----------------------------------------------------------------------------
  if (boot == FALSE) {
    # Estimate of standard error
    se.dr.att <- stats::sd(dr.att.inf.func)/sqrt(n)
    # Estimate of upper boudary of 95% CI
    uci <- dr.att + qnorm(1-(1-alpha)/2) * se.dr.att
    # Estimate of lower doundary of 95% CI
    lci <- dr.att - qnorm(1-(1-alpha)/2) * se.dr.att
    #Create this null vector so we can export the bootstrap draws too.
    dr.boot <- NULL
  }
  
  if (boot == TRUE) {
    if (is.null(nboot) == TRUE) nboot = 999
    # default to "multiplier bootstrap" 
    # do multiplier bootstrap
    dr.boot <- mboot.did(dr.att.inf.func, nboot)
    # get bootstrap std errors based on IQR
    se.dr.att <- stats::IQR(dr.boot) / (stats::qnorm(0.75) - stats::qnorm(0.25))
    # get symmtric critival values
    cv <- stats::quantile(abs(dr.boot/se.dr.att), probs = alpha)
    # Estimate of upper boudary of 95% CI
    uci <- dr.att + cv * se.dr.att
    # Estimate of lower doundary of 95% CI
    lci <- dr.att - cv * se.dr.att
    
  }
  
  
  if(inffunc == FALSE) dr.att.inf.func <- NULL
  #---------------------------------------------------------------------
  # record the call
  call.param <- match.call()
  # Record all arguments used in the function
  argu <- mget(names(formals()), sys.frame(sys.nframe()))
  boot <- ifelse(argu$boot == TRUE, TRUE, FALSE)
  argu <- list(
    stabilized = stabilized,
    boot = boot,
    #boot.type = boot.type,
    nboot = nboot,
    type = "dr,nonstnr"
  )
  ret <- (list(ATT = dr.att,
               se = se.dr.att,
               uci = uci,
               lci = lci,
               boots = dr.boot,
               att.inf.func = dr.att.inf.func,
               call.param = call.param,
               argu = argu))
  
  # Define a new class
  class(ret) <- "drdid"
  # return the list
  return(ret)
}



#' Stationary DR DiD estimator (Sant'Anna & Zhao, 2020)
#'
#' Uses two-group propensity score (treated vs control) under the assumption
#' that covariate distributions are stationary across periods.
#'
#' @inheritParams drdid_nonstnr
#' @return List with ATT, se, uci, lci, boots, att.inf.func.
drdid_stnr <- function(y, d, post, fit.ps, fit.or,
                       stabilized = TRUE, i.weights = NULL,
                       boot = FALSE, nboot = NULL,
                       inffunc = FALSE, alpha = 0.95) {
  # D as vector
  d <- as.vector(d)
  # Sample size
  n <- length(d)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  
  # Get the PS estimates
  ps.fit.treat = fit.ps[,1] + fit.ps[,2]
  
  pmin <- 0.005
  ps.fit.treat[ps.fit.treat < pmin] = pmin
  ps.fit.treat[ps.fit.treat > (1-pmin)] = 1-pmin
  
  
  # Get the OR estimates
  out.y.treat.post = fit.or[,1]
  out.y.treat.pre  = fit.or[,2]
  out.y.cont.post  = fit.or[,3]
  out.y.cont.pre   = fit.or[,4]
  
  #-----------------------------------------------------------------------------
  # Compute the ATT
  # First, the weights
  w.treat.post <-  d * post
  w.treat.pre  <-  d * (1 - post) 
  w.cont.post  <-  (1 - d) * post * ps.fit.treat / (1-ps.fit.treat)
  w.cont.pre   <-  (1 - d) * (1 - post) * ps.fit.treat / (1-ps.fit.treat)
  
  # Elements of the influence function (summands)
  if (stabilized == TRUE)
  {
    eta.treat.post <- w.treat.post * (y - out.y.treat.post) / mean(w.treat.post)
    eta.treat.pre  <- w.treat.pre * (y - out.y.treat.pre) / mean(w.treat.pre)
    eta.cont.post  <- w.cont.post * (y - out.y.cont.post) / mean(w.cont.post)
    eta.cont.pre   <- w.cont.pre  * (y - out.y.cont.pre) / mean(w.cont.pre)
    eta.tau.x  <- d * (out.y.treat.post -  out.y.treat.pre -
                         out.y.cont.post  +  out.y.cont.pre) / mean(d)
  } else 
  {
    eta.treat.post <- w.treat.post * (y - out.y.treat.post) / mean(w.treat.post)
    eta.treat.pre  <- w.treat.pre * (y - out.y.treat.pre) / mean(w.treat.pre)
    eta.cont.post  <- w.cont.post * (y - out.y.cont.post) / mean(w.treat.post)
    eta.cont.pre   <- w.cont.pre  * (y - out.y.cont.pre) / mean(w.treat.pre)
    eta.tau.x  <- d * (out.y.treat.post -  out.y.treat.pre -
                         out.y.cont.post  +  out.y.cont.pre) / mean(d)
  }
  
  
  
  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)
  att.tau.x <-  mean(eta.tau.x)
  # ATT estimator
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre)  + att.tau.x
  
  #-----------------------------------------------------------------------------
  # get the influence function and compute standard error
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # Influence function for the treated component
  inf.treat <- eta.treat.post - eta.treat.pre
  #-----------------------------------------------------------------------------
  # Influence function for the control component
  inf.cont <- eta.cont.post - eta.cont.pre
  #-----------------------------------------------------------------------------
  #get the influence function of the DR estimator (put all pieces together)
  dr.att.inf.func1 <- inf.treat - inf.cont
  #-----------------------------------------------------------------------------
  # Now, we only need to get the influence function of the adjustment terms
  # First, the terms as if all OR parameters were known
  inf.eff <- eta.tau.x - d * dr.att / mean(d)
  #-----------------------------------------------------------------------------
  #get the influence function of the locally efficient DR estimator (put all pieces together)
  dr.att.inf.func <- dr.att.inf.func1 + inf.eff
  #-----------------------------------------------------------------------------
  if (boot == FALSE) {
    # Estimate of standard error
    se.dr.att <- stats::sd(dr.att.inf.func)/sqrt(n)
    # Estimate of upper boudary of 95% CI
    uci <- dr.att + qnorm(1-(1-alpha)/2) * se.dr.att
    # Estimate of lower doundary of 95% CI
    lci <- dr.att - qnorm(1-(1-alpha)/2) * se.dr.att
    #Create this null vector so we can export the bootstrap draws too.
    dr.boot <- NULL
  }
  
  if (boot == TRUE) {
    if (is.null(nboot) == TRUE) nboot = 999
    # default to "multiplier bootstrap" 
    dr.boot <- mboot.did(dr.att.inf.func, nboot)
    # get bootstrap std errors based on IQR
    se.dr.att <- stats::IQR(dr.boot) / (stats::qnorm(0.75) - stats::qnorm(0.25))
    # get symmtric critival values
    cv <- stats::quantile(abs(dr.boot/se.dr.att), probs = alpha)
    # Estimate of upper boudary of 95% CI
    uci <- dr.att + cv * se.dr.att
    # Estimate of lower doundary of 95% CI
    lci <- dr.att - cv * se.dr.att
    
  }
  
  if(inffunc == FALSE) dr.att.inf.func <- NULL
  #---------------------------------------------------------------------
  # record the call
  call.param <- match.call()
  # Record all arguments used in the function
  argu <- mget(names(formals()), sys.frame(sys.nframe()))
  boot <- ifelse(argu$boot == TRUE, TRUE, FALSE)
  argu <- list(
    stabilized = stabilized,
    boot = boot,
    nboot = nboot,
    type = "dr,stnr"
  )
  ret <- (list(ATT = dr.att,
               se = se.dr.att,
               uci = uci,
               lci = lci,
               boots = dr.boot,
               att.inf.func = dr.att.inf.func,
               call.param = call.param,
               argu = argu))
  
  # Define a new class
  class(ret) <- "drdid"
  # return the list
  return(ret)
}



#' TWFE DiD estimator for repeated cross-sections
#'
#' Estimates ATT via OLS: y ~ d*post + post + d + covariates.
#' The coefficient on d*post is the TWFE ATT estimate.
#'
#' @param y          Outcome vector.
#' @param post       Post-period indicator.
#' @param d          Treatment indicator.
#' @param covariates Optional covariate matrix.
#' @param i.weights  Optional observation weights.
#' @param boot       Logical; use bootstrap.
#' @param boot.type  Bootstrap type: "weighted" or "multiplier".
#' @param nboot      Number of bootstrap draws.
#' @param inffunc    Logical; return influence function.
#' @param alpha      Nominal coverage level(s).
#' @return List with ATT, se, uci, lci, boots, att.inf.func.
twfe_did_rc <- function(y, post, d, covariates = NULL, i.weights = NULL,
                        boot = FALSE, boot.type = "weighted", nboot = NULL,
                        inffunc = FALSE, alpha = 0.95) {
  #-----------------------------------------------------------------------------
  # D as vector
  d <- as.vector(d)
  # Sample size
  n <- length(d)
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  #-----------------------------------------------------------------------------
  #Create dataset for TWFE approach
  if (is.null(covariates)) {
    x = NULL
  } else {
    if(all(as.matrix(covariates)[,1] == rep(1,n))) {
      # Remove intercept if you include it
      covariates <- as.matrix(covariates)
      covariates <- covariates[,-1]
      if(dim(covariates)[2]==0) {
        covariates = NULL
        x = NULL
      }
    }
  }
  
  if(!is.null(covariates))  x <- as.matrix(covariates)
  dd <- d
  post <- post
  i.weights <- as.vector(i.weights)
  #---------------------------------------------------------------------------
  #Estimate TWFE regression
  if(!is.null(x)){
    reg <- stats::lm(y ~  dd:post + post + dd + x, x = TRUE, weights = i.weights)
  }
  if(is.null(x)){
    reg <- stats::lm(y ~  dd:post + post + dd, x = TRUE, weights = i.weights)
  }
  twfe.att <- reg$coefficients["dd:post"]
  #-----------------------------------------------------------------------------
  #Elemenets for influence functions
  inf.reg <- (i.weights * reg$x * reg$residuals) %*%
    base::solve(crossprod(i.weights * reg$x, reg$x) / dim(reg$x)[1])
  
  sel.theta <- matrix(c(rep(0, dim(inf.reg)[2])))
  
  index.theta <- which(dimnames(reg$x)[[2]]=="dd:post",
                       arr.ind = TRUE)
  
  sel.theta[index.theta, ] <- 1
  #-----------------------------------------------------------------------------
  #get the influence function of the TWFE regression
  twfe.inf.func <- as.vector(inf.reg %*% sel.theta)
  #-----------------------------------------------------------------------------
  if (boot == FALSE) {
    # Estimate of standard error
    se.twfe.att <- stats::sd(twfe.inf.func)/sqrt(length(twfe.inf.func))
    # Estimate of upper boudary of 95% CI
    uci <- twfe.att + qnorm(1-(1-alpha)/2) * se.twfe.att
    # Estimate of lower doundary of 95% CI
    lci <- twfe.att - qnorm(1-(1-alpha)/2) * se.twfe.att
    #Create this null vector so we can export the bootstrap draws too.
    twfe.boot <- NULL
  }
  
  if (boot == TRUE) {
    if (is.null(nboot) == TRUE) nboot = 999
    if(boot.type == "multiplier"){
      # do multiplier bootstrap
      twfe.boot <- mboot.did(twfe.inf.func, nboot)
      # get bootstrap std errors based on IQR
      se.twfe.att <- stats::IQR(twfe.boot) / (stats::qnorm(0.75) - stats::qnorm(0.25))
      # get symmtric critival values
      cv <- stats::quantile(abs(twfe.boot/se.twfe.att), probs = alpha)
      # Estimate of upper boudary of 95% CI
      uci <- twfe.att + cv * se.twfe.att
      # Estimate of lower doundary of 95% CI
      lci <- twfe.att - cv * se.twfe.att
    } else {
      # do weighted bootstrap
      twfe.boot <- unlist(lapply(1:nboot, wboot_twfe_rc,
                                 n = n, y = y, dd = dd, post = post, x = x, i.weights = i.weights))
      # get bootstrap std errors based on IQR
      se.twfe.att <- stats::IQR((twfe.boot - twfe.att)) / (stats::qnorm(0.75) - stats::qnorm(0.25))
      # get symmtric critival values
      cv <- stats::quantile(abs((twfe.boot - twfe.att)/se.twfe.att), probs = alpha)
      # Estimate of upper boudary of 95% CI
      uci <- twfe.att + cv * se.twfe.att
      # Estimate of lower doundary of 95% CI
      lci <- twfe.att - cv * se.twfe.att
      
    }
  }
  
  
  if(inffunc == FALSE) att.inf.func <- NULL
  return(list(ATT = twfe.att,
              se = se.twfe.att,
              uci = uci,
              lci = lci,
              boots = twfe.boot,
              att.inf.func = twfe.inf.func))
}



#' Mammen (1993) multiplier bootstrap for DiD influence functions
#'
#' @param linrep Influence function vector (length n).
#' @param nboot  Number of bootstrap draws.
#' @return Vector of bootstrapped ATT deviations.
mboot.did <- function(linrep, nboot) {
  # Use the Mammen (1993) binary V's
  k1 = 0.5 * (1 - 5^0.5)
  k2 = 0.5 * (1 + 5^0.5)
  pkappa = 0.5 * (1 + 5^0.5)/(5^0.5)
  
  n <- length(linrep)
  bootapply <- function(nn,  n = n, linrep = linrep) {
    v <- stats::rbinom(n, 1, pkappa)
    v <- ifelse(v == 1, k1, k2)
    b.did <- mean(linrep * v)
    return(b.did)
  }
  
  boot.did <- unlist(lapply(1:nboot, bootapply, n = n, linrep = linrep))
  
  
}


#' Cluster-level multiplier bootstrap for a single DiD estimator
#'
#' Draws Mammen weights at the cluster (group) level.
#'
#' @param linrep Influence function vector (length n).
#' @param group  Cluster group identifiers (length n).
#' @param nboot  Number of bootstrap draws.
#' @return Vector of bootstrapped ATT deviations.
boot.did.cluster <- function(linrep, group, nboot) {
  # Mammen (1993) binary V weights at cluster level
  k1 <- 0.5 * (1 - 5^0.5)
  k2 <- 0.5 * (1 + 5^0.5)
  pkappa <- 0.5 * (1 + 5^0.5) / (5^0.5)
  
  
  bootgroup <- data.frame(group)
  
  bootapply <- function(nn) {
    # seed.run <- Seed[nn,]
    # set.seed(seed.run, "L'Ecuyer-CMRG") ## to make each run fully reproducible
    bootgroup <- bootgroup %>%
      group_by(group, .add=T) %>%
      mutate(v = stats::rbinom(1, 1, pkappa)) %>%
      ungroup()

    v <- ifelse(bootgroup$v == 1, k1, k2)
    b.did <- mean(linrep * v)
  }
  
  boot.did <- unlist(lapply(1:nboot, bootapply))
  
}


#' Cluster-level multiplier bootstrap for joint DR DiD (stnr + nstnr)
#'
#' Jointly bootstraps both the stationary and non-stationary estimators,
#' and computes bootstrapped Wald statistics for the Hausman test.
#'
#' @param linrep n x 2 matrix of influence functions (col 1: stnr, col 2: nstnr).
#' @param group  Cluster group identifiers (length n).
#' @param nboot  Number of bootstrap draws.
#' @return List with b.stnr, b.nstnr (bootstrap ATT draws), wald.stat.
boot.drdid.cluster <- function(linrep, group, nboot) {
  # Mammen (1993) binary V weights at cluster level
  n <- length(group)
  k1 <- 0.5 * (1 - 5^0.5)
  k2 <- 0.5 * (1 + 5^0.5)
  pkappa <- 0.5 * (1 + 5^0.5) / (5^0.5)
  
  bootgroup <- data.frame(group)
  
  bootapply <- function(nn) {
    bootgroup <- bootgroup %>%
      group_by(group, .add=T) %>%
      mutate(v = stats::rbinom(1, 1, pkappa)) %>%
      ungroup()
    v <- ifelse(bootgroup$v == 1, k1, k2)
    
    return(list(se = colMeans(linrep * v),
                wald  = n* (mean(linrep[,1]*v - linrep[,2]*v))^2/mean((linrep[,1]*v-linrep[,2]*v)^2)   
    )
    )
  }
  # first col of linrep is stnr, second is nstnr
  boot.did <- matrix(unlist(lapply(1:nboot, bootapply)),  ncol = 3, byrow = TRUE)
  
  return(list(b.stnr = boot.did[,1],
              b.nstnr = boot.did[,2],
              wald.stat  = boot.did[,3]
  )
  )
  
}
