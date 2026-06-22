#' Non-stationary doubly robust DiD estimator (repeated cross-sections)
#'
#' Estimates the ATT allowing the covariate distribution to change across
#' time/groups (Sant'Anna & Xu), together with its influence function,
#' standard error, and confidence interval.
#'
#' @param y Outcome vector.
#' @param d Treated-group indicator (1 = treated group).
#' @param post Post-period indicator (1 = post period).
#' @param fit.ps Matrix of fitted generalized propensity scores, columns
#'   (treat-post, treat-pre, control-post, control-pre).
#' @param fit.or Matrix of fitted outcome regressions, same four columns.
#' @param stabilized Logical; use stabilized weights (default TRUE).
#' @param i.weights Optional non-negative sampling weights.
#' @param boot Logical; if TRUE use a multiplier bootstrap for inference.
#' @param nboot Number of bootstrap draws (default 999 when boot = TRUE).
#' @param inffunc Logical; if TRUE return the influence function.
#' @return Object of class "drdid": ATT, se, uci, lci, boots,
#'   att.inf.func, call.param, argu.
drdid_nonstnr = function(y, d, post, fit.ps, fit.or, 
                         stabilized = TRUE, i.weights = NULL, 
                         boot = FALSE, nboot = NULL,
                         inffunc = TRUE){
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
  
  # Get the OR estimates
  ps.fit.treat.post = fit.ps[,1]
  ps.fit.treat.pre  = fit.ps[,2]
  ps.fit.cont.post  = fit.ps[,3]
  ps.fit.cont.pre   = fit.ps[,4]
  ps.fit.treat = ps.fit.treat.post + ps.fit.treat.pre
  
  # Get the PS estimates
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
    uci <- dr.att + 1.96 * se.dr.att
    # Estimate of lower doundary of 95% CI
    lci <- dr.att - 1.96 * se.dr.att
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
    cv <- stats::quantile(abs(dr.boot/se.dr.att), probs = 0.95)
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
               #ps.flag = pscore.ipt$flag,
               att.inf.func = dr.att.inf.func,
               call.param = call.param,
               argu = argu))
  
  # Define a new class
  class(ret) <- "drdid"
  # return the list
  return(ret)
}



#' Stationary doubly robust DiD estimator (repeated cross-sections)
#'
#' Estimates the ATT under a stationary covariate distribution
#' (Sant'Anna & Zhao, 2020), with influence function, standard error, and
#' confidence interval. Arguments and return value mirror drdid_nonstnr().
#'
#' @param y Outcome vector.
#' @param d Treated-group indicator (1 = treated group).
#' @param post Post-period indicator (1 = post period).
#' @param fit.ps Matrix of fitted generalized propensity scores, columns
#'   (treat-post, treat-pre, control-post, control-pre).
#' @param fit.or Matrix of fitted outcome regressions, same four columns.
#' @param stabilized Logical; use stabilized weights (default TRUE).
#' @param i.weights Optional non-negative sampling weights.
#' @param boot Logical; if TRUE use a multiplier bootstrap for inference.
#' @param nboot Number of bootstrap draws (default 999 when boot = TRUE).
#' @param inffunc Logical; if TRUE return the influence function.
#' @return Object of class "drdid": ATT, se, uci, lci, boots,
#'   att.inf.func, call.param, argu.
drdid_stnr = function(y, d, post, fit.ps, fit.or, 
                      stabilized = TRUE, i.weights = NULL, 
                      boot = FALSE, nboot = NULL,
                      inffunc = TRUE){
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
  
  # Get the OR estimates
  ps.fit.treat.post = fit.ps[,1]
  ps.fit.treat.pre  = fit.ps[,2]
  ps.fit.cont.post  = fit.ps[,3]
  ps.fit.cont.pre   = fit.ps[,4]
  ps.fit.treat = ps.fit.treat.post + ps.fit.treat.pre
  
  # Get the PS estimates
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
  dr.att.inf.func <- inf.eff  + dr.att.inf.func1 
  #-----------------------------------------------------------------------------
  if (boot == FALSE) {
    # Estimate of standard error
    se.dr.att <- stats::sd(dr.att.inf.func)/sqrt(n)
    # Estimate of upper boudary of 95% CI
    uci <- dr.att + 1.96 * se.dr.att
    # Estimate of lower doundary of 95% CI
    lci <- dr.att - 1.96 * se.dr.att
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
    cv <- stats::quantile(abs(dr.boot/se.dr.att), probs = 0.95)
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


#' Two-way fixed-effects DiD estimator (repeated cross-sections)
#'
#' OLS DiD with a treatment-by-post interaction (the TWFE benchmark), with an
#' influence-function-based standard error and confidence interval.
#'
#' @param y Outcome vector.
#' @param post Post-period indicator (1 = post period).
#' @param d Treated-group indicator (1 = treated group).
#' @param covariates Optional covariate design matrix (a leading intercept
#'   column, if present, is dropped).
#' @param i.weights Optional non-negative sampling weights.
#' @param boot Logical; if TRUE bootstrap for inference.
#' @param boot.type "weighted" (exponential multiplier) or "multiplier" (Mammen).
#' @param nboot Number of bootstrap draws (default 999 when boot = TRUE).
#' @param inffunc Logical; if TRUE return the influence function.
#' @return List: ATT, se, uci, lci, boots, att.inf.func.
twfe_did_rc <- function(y, post, d, covariates = NULL, i.weights = NULL,
                        boot = FALSE, boot.type = "weighted", nboot = NULL,
                        inffunc = FALSE){
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
    uci <- twfe.att + 1.96 * se.twfe.att
    # Estimate of lower doundary of 95% CI
    lci <- twfe.att - 1.96 * se.twfe.att
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
      cv <- stats::quantile(abs(twfe.boot/se.twfe.att), probs = 0.95)
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
      cv <- stats::quantile(abs((twfe.boot - twfe.att)/se.twfe.att), probs = 0.95)
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
              att.inf.func = att.inf.func))
}


# Bootrstapped Doubly Robust Difference-in-Differences
# 2 periods and 2 groups
#' Multiplier (Mammen) bootstrap of a DiD influence function
#'
#' @param linrep Influence-function (linear representation) vector.
#' @param nboot Number of bootstrap draws.
#' @return Vector of nboot bootstrapped DiD estimates.
mboot.did = function(linrep, nboot){
  # Use the Mammen(1993) binary V's
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



#' Single weighted-bootstrap replication of the TWFE estimator
#'
#' One exponential-multiplier bootstrap draw of twfe_did_rc(), for use inside
#' lapply().
#'
#' @param nn Replication index (ignored; required by lapply()).
#' @param n Sample size.
#' @param y Outcome vector.
#' @param dd Treated-group indicator.
#' @param post Post-period indicator.
#' @param x Covariate matrix (or NULL).
#' @param i.weights Sampling weights.
#' @return Bootstrapped TWFE ATT for a single draw.
wboot_twfe_rc <- function(nn, n, y, dd, post, x, i.weights){
  #-----------------------------------------------------------------------------
  v <- stats::rexp(n)
  #v <- v / mean(v)
  #weights for the bootstrap
  b.weights <- as.vector(i.weights * v)
  #Compute the TWFE Regression
  if(!is.null(x)){
    reg.b <- stats::lm(y ~  dd:post + post + dd + x, weights = b.weights)
  } else{
    reg.b <- stats::lm(y ~  dd:post + post + dd, weights = b.weights)
  }
  twfe.att.b <- reg.b$coefficients["dd:post"]
  #-----------------------------------------------------------------------------
  return(twfe.att.b)
}
