#' Data-generating process for DRDID simulations
#'
#' @param delta Mixture weight controlling deviation from the NCC condition.
#' @param n Sample size.
#' @param xi_ps Scaling factor for propensity score coefficients.
#' @return List with simulated outcomes, covariates, and ground-truth targets.
dgps_did <- function(delta, n, xi_ps = 1) {
  # Generate covariates
  x1 <- 2*(stats::runif(n) - 0.5)
  x2 <- 2*(stats::runif(n) - 0.5)
  x3 <- stats::rbinom(n, 1, prob = 0.5)
  x4 <- stats::rbinom(n, 1, prob = 0.5)
  x5 <- stats::rbinom(n, 3, prob = 0.5)
  x6 <- stats::rbinom(n, 3, prob = 0.5)
  x  <- cbind(x1, x2, x3, x4, x5, x6)
  dim_x <- c(2,2,2)
  
  # Generate treatment groups  
  gamma_11 <- c(0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0)
  
  gamma_10 <- xi_ps * c(0, 0.4, 0.4, -0.4, -0.4, 0,
                        0.2, 0.2, 0.2, 0.2, 0.1, -0.1,
                        0.1, 0.1, 0.1, 0.1,
                        -0.1, -0.1, -0.1, -0.1,
                        0.1, -0.1, -0.1, 0.1)
  
  gamma_01 <- xi_ps * c(0, 0.8, 0.4, 0.4, -0.4, 0.4,
                        0.2, -0.2, 0.2, -0.2, 0, 0,
                        0, 0, 0, 0,
                        0.1, 0.1, 0.1, 0.1,
                        0.0, 0.1, 0.0, 0.1)
  
  gamma_00 <- xi_ps * c(0, 0.4, 0.8, -0.4, 0.4, -0.4,
                        -0.2, 0.2, -0.2, 0.2, 0, 0,
                        0.1, 0.1, 0.1, 0.1,
                        0, 0, 0, 0,
                        0.1, 0.0, 0.1, 0)
  
  
  # Generalized Propensity Scores
  design.x <- as.matrix(cbind(1, x1, x2, x1^2, x2^2, x1*x2, 
                              x3, x4, x5, x6, x3*x4, x5*x6,
                              x1*x3, x1*x4, x1*x5, x1*x6,
                              x2*x3, x2*x4, x2*x5, x2*x6,
                              x3*x5, x3*x6, x4*x5, x4*x6))
  ps.11 <- exp(design.x %*% gamma_11)
  ps.10 <- exp(design.x %*% gamma_10)
  ps.01 <- exp(design.x %*% gamma_01)
  ps.00 <- exp(design.x %*% gamma_00)
  ps.denom <- ps.11 + ps.10 + ps.01 + ps.00 
  
  
  if (delta >= 0 && delta <= 1) { 
    # Power analysis: mixture of DGPs 1 and 2 with weights governed by delta
    pr.treat <- (ps.10 + ps.11) / ps.denom
    pr.post <- 0.467
    
    p10 <- (1 - pr.post) * pr.treat * delta  +  ps.10 / ps.denom * (1 - delta)
    p01 <- pr.post * (1 - pr.treat) * delta  +  ps.01 / ps.denom * (1 - delta)
    p00 <- (1 - pr.post) * (1 - pr.treat) * delta  +  ps.00 / ps.denom * (1 - delta)
    p11 <- 1 - p10 - p01 - p00
    
    pp10 = (1-pr.post) * pr.treat 
    pp01 = pr.post * (1-pr.treat) 
    pp00 = (1-pr.post) * (1-pr.treat) 
    pp11 = 1 - pp10 - pp01 - pp00
    
    u  <- runif(n)
    DT <- (u <= p10) + 2 * (u > p10 & u <= (p10 + p01)) + 3 * (u > (p10 + p01) & u <= (p10 + p01 + p00))
    d  <- as.numeric((DT == 0) + (DT  == 1))
    post <- as.numeric((DT == 0) + (DT == 2))
    DT_mat <- cbind(d*post, d*(1-post), (1-d)*post, (1-d)*(1-post))
    
    # Generate aux indexes for the potential outcomes
    baseline  <- 27.4*x1 + 27.4*x2  + 13.7*x1^2  + 13.7*x1*x2 + 13.7*x2^2
    index.lin <- 210 + baseline
    # gen treatment heterogeneity
    v1 <- baseline
    v0 <- 0
    index.unobs.het <- d * v1 + (1-d) * v0
    eps <- stats::rnorm(n, mean = 0, sd = 1)
    v <-  index.unobs.het + eps
    
    # gen time trend
    index.trend <-  baseline
    # gen treatment effect
    index.att <-  27.4*x1 + 13.7*x2 + 13.7*(x5+x3+x4+x6)/2 -15
    
    # gen realized outcome at time 0
    y10 <- index.lin + v + stats::rnorm(n)
    y00 <- index.lin + v + stats::rnorm(n)
    # gen outcomes at time 1
    # first let's generate potential outcomes: y_1_potential
    y01 <- index.lin + v + stats::rnorm(n) + #This is the baseline
      index.trend #this is for the trend based on X
    y11 <- index.lin + v + stats::rnorm(n) + #This is the baseline
      index.trend + #this is for the trend based on X
      index.att # This is the treatment effects
    
    # Get the conditional mean functions
    m10 <- index.lin + v1
    m00 <- index.lin + v0
    m01 <- index.lin + index.trend + v0
    m11 <- index.lin + index.trend + index.att + v1
    
    # Get infeasible att
    att.unf <- (mean(d *post * y11) - mean(d *post * y01) -  (mean(d * post * y10) - mean(d *post * y00)))/mean(post * d)
    att.unf.stnr <- (mean(d * y11) - mean(d * y01) -  (mean(d * y10) - mean(d * y00)))/mean( d)
    
    # Set true parameter values and semiparametric efficiency bound
    if (delta == 0){
      att.true <- 4.308165
      eff <- 1753.598
    } else if (delta == 1){
      att.true <- 9.128423
      eff <- 796.8117
    } else if (delta == 0.5){
      att.true <- 6.819859
      eff <- 1731.674
    } else{
      att.true <- NA_real_
      eff <- NA_real_
    }
    
    # Generate realized outcome at time 0 and 1
    y.post <- d * y11 + (1 - d) * y01
    y.pre <- d * y10 + (1 - d) * y00
    # Generate observed outcome
    y <- post * y.post + (1 - post) * y.pre
    
  } else {
    stop(sprintf("Invalid delta (%s): must lie in the interval [0, 1].", delta))
  }
  
  xcov <-   as.matrix(cbind(x1, x2, x1^2, x2^2, x1*x2, 
                            x3, x4, x5, x6, x3*x4, x5*x6,
                            x1*x3, x1*x4, x1*x5, x1*x6,
                            x2*x3, x2*x4, x2*x5, x2*x6,
                            x3*x5, x3*x6, x4*x5, x4*x6))
  
  return(list(y10 = y10, y11 = y11, y01 = y01, y00 = y00,
              m10 = m10, m11 = m11, m01 = m01, m00 = m00,
              m = cbind(m11, m10, m01, m00),
              y = y, d = d, 
              post = post,
              dpost = DT_mat,
              dty   = DT_mat * y,
              covariates = x, 
              covariates_sat  = xcov,
              dim_covariates = dim_x,
              pi = cbind(p11, p10, p01, p00),
              pp = cbind(pp11, pp10, pp01, pp00),
              att.true = att.true, 
              att.unf = att.unf,
              att.unf.stnr = att.unf.stnr,
              eff = eff))
  
}



#' True ATT and semiparametric efficiency bound for a given DGP
#'
#' Computes, by large-sample Monte Carlo, the true ATT and the semiparametric
#' efficiency bound implied by the population influence function, for either
#' the non-stationary or the stationary case.
#'
#' @param delta Mixture weight controlling deviation from the NCC condition.
#' @param n Sample size for the Monte Carlo approximation (use large n).
#' @param stnr Logical; if TRUE use the stationary efficiency bound.
#' @param Xsi.ps Scaling factor for propensity score coefficients.
#' @return List with att_true (true ATT) and eff (efficiency bound).
se.bound <- function(delta, n, stnr = FALSE, Xsi.ps= 1){
  data  = dgps_did(delta, n, xi_ps = 1)
  
  p11 = data$pi[,1]
  p10 = data$pi[,2]
  p01 = data$pi[,3]
  p00 = data$pi[,4]
  pi = p11 + p10
  d    = data$d
  post = data$post
  y    = data$y
  
  m11 =   data$m11 
  m10 =   data$m10 
  m01 =   data$m01 
  m00 =   data$m00 
  
  if (stnr == FALSE){
    eta.treat.post <- d * post * (y - m11)/mean(d * post)
    eta.treat.pre  <- d * (1-post) * p11 * (y - m10)/(p10 * mean(d * post ))
    eta.cont.post  <- (1-d) * post * p11 * (y - m01)/(p01 * mean(d * post ))
    eta.cont.pre   <- (1-d) * (1-post) * p11 * (y - m00)/(p00 * mean(d * post ))
    # Influence function for the treated component
    inf.treat <- eta.treat.post - eta.treat.pre
    #-----------------------------------------------------------------------------
    # Influence function for the control component
    inf.cont <- eta.cont.post - eta.cont.pre
    #-----------------------------------------------------------------------------
    #get the influence function of the DR estimator (put all pieces together)
    dr.att.inf.func1 <- inf.treat - inf.cont
    #-----------------------------------------------------------------------------
    eta.tau.x  <- d * post * (m11-m10-m01+m00) / mean(d * post)
    # Estimator of each component
    att.treat.pre <- mean(eta.treat.pre)
    att.treat.post <- mean(eta.treat.post)
    att.cont.pre <- mean(eta.cont.pre)
    att.cont.post <- mean(eta.cont.post)
    att.tau.x <-  mean(eta.tau.x)
    # estimate ATT 
    att_true <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre)  + att.tau.x
    att_unf <- (mean(d *post * data$y11) - mean(d *post * data$y01) -  (mean(d * post * data$y10) - mean(d *post * data$y00)))/mean(post * d)
    #att_stnr <- (mean(d * data$y11) - mean(d * data$y01) -  (mean(d *  data$y10) - mean(d *data$y00)))/mean(d)
    # Now, we only need to get the influence function of the adjustment terms
    # First, the terms as if all OR parameters were known
    inf.eff <- eta.tau.x - d * post * att_true / mean(d * post)
  } else{
    eta.treat.post <- d * post * (y - m11)/mean(d * post)
    eta.treat.pre  <- d * (1-post) * p11 * (y - m10)/(p10 * mean(d * post ))
    eta.cont.post  <- (1-d) * post * pi * (y - m01)/(pi * mean(d * post ))
    eta.cont.pre   <- (1-d) * (1-post) * pi * (y - m00)/((1-pi) * mean(d * post ))
    # Influence function for the treated component
    inf.treat <- eta.treat.post - eta.treat.pre
    #-----------------------------------------------------------------------------
    # Influence function for the control component
    inf.cont <- eta.cont.post - eta.cont.pre
    #-----------------------------------------------------------------------------
    #get the influence function of the DR estimator (put all pieces together)
    dr.att.inf.func1 <- inf.treat - inf.cont
    #-----------------------------------------------------------------------------
    eta.tau.x  <- d *  (m11-m10-m01+m00) / mean(d)
    
    # Estimator of each component
    att.treat.pre <- mean(eta.treat.pre)
    att.treat.post <- mean(eta.treat.post)
    att.cont.pre <- mean(eta.cont.pre)
    att.cont.post <- mean(eta.cont.post)
    att.tau.x <-  mean(eta.tau.x)
    # estimate ATT 
    att_true <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre)  + att.tau.x
    att_unf <- (mean(d * data$y11) - mean(d * data$y01) -  (mean(d *  data$y10) - mean(d *data$y00)))/mean(d)
    # Now, we only need to get the influence function of the adjustment terms
    # First, the terms as if all OR parameters were known
    inf.eff <- eta.tau.x - d * att_true / mean(d)
  }
  
  #-----------------------------------------------------------------------------
  #get the influence function of the locally efficient DR estimator (put all pieces together)
  dr.att.inf.func <-  inf.eff + dr.att.inf.func1
  
  eff <- stats::var(dr.att.inf.func)
  
  return(list(att_true = att_true,
              eff = eff))
}