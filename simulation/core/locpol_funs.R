# Constants related to initial values
varc_init = 0.7
vard_init = 0.3
lbc_init = 0.1
ubc_init = 2.0
lbd_init = 0.1
ubd_init = 0.9

#-----------------------------------------------------------------------------
# Cross-validated Bandwidth for the GPS model
#' Cross-validated bandwidth for the GPS local-polynomial model
#'
#' @param dp Data list with dpost, covariates, and dim_covariates.
#' @param cv_method CV criterion ("cv.ml" likelihood-based or "cv.ls"
#'   least-squares).
#' @param lp_type Local-polynomial link (e.g. "logit").
#' @param n_start Number of random restarts for the bandwidth optimizer.
#' @return Cross-validation bandwidth object returned by cv_ps().
locpol.ps.bw <- function(dp,
                         cv_method, 
                         lp_type,
                         n_start){
 

  #---------------------------------------------------------------------  
  # initialize the bandwidth parameters 
  if (n_start > 0){
    
    rc_bw = runif(n_start-1)
    bwc_init =  c(varc_init, (ubc_init - lbc_init)*rc_bw+lbc_init)
    ru_bw = runif(n_start-1)
    bwu_init =  c(vard_init, (ubd_init - lbd_init)*ru_bw+lbd_init)
    ro_bw = runif(n_start-1)
    bwo_init =  c(vard_init, (ubd_init - lbd_init)*ro_bw+lbd_init)
    
    # put them together
    bw_init = rbind(bwc_init, bwu_init, bwo_init)
    
    } else {
      stop("number of restarts must be a positive integer.")
    }
  
  #---------------------------------------------------------------------  
  # estimate the cross-validated bandwidth for nonparametric PS model
  cv.bw <- cv_ps(dp$dpost, 
                 dp$covariates, 
                 dp$dim_covariates, 
                 cv_method, 
                 lp_type, 
                 n_start,
                 bw_init)
  
  return(cv.bw)
    
}

#-----------------------------------------------------------------------------
# GPS estimates
#' Estimate the generalized propensity score by local-polynomial logit
#'
#' Selects the bandwidth by cross-validation when none is supplied, then fits
#' the GPS and trims fitted values away from 0 and 1.
#'
#' @param dp Data list with dpost, covariates, and dim_covariates.
#' @param bws Optional bandwidth vector; if NULL, chosen by CV.
#' @param cv_method CV criterion ("cv.ml" or "cv.ls").
#' @param lp_type Local-polynomial link (e.g. "logit").
#' @param list_control List with n_start (restarts) and ps_min (trimming bound).
#' @return List: fitted.values, bws, call.param, argu.
locpol.ps.fit <- function(dp,
                          bws = NULL, 
                          cv_method= "cv.ml", 
                          lp_type= "logit",
                          list_control = list(n_start = 3,
                                              ps_min = 1e-5)){
  
  #---------------------------------------------------------------------  
  # check if bandwidth is provided, if not, estimate it by cross validation
  if(is.null(bws)) {
     bws_list <- locpol.ps.bw(dp,
                              cv_method = cv_method,
                              lp_type   = lp_type,
                              n_start   = list_control$n_start
                             )
    
      bws = as.numeric(bws_list$bw_cv)
  } 
  
  # # check if bws has the correct dimension
  #  if ( !all(which(bws != 0) == which(dp$dim_covariates != 0)) ) {
  #    stop("must provide at least one bandwidth for each type of covariates in the data.")
  #  }
  
  #---------------------------------------------------------------------
  # compute the generalized propensity scores by logistic local polynomial regression

  pscore.fitted   <-  locpolfit_ps(bws, 
                                   dp$dpost, 
                                   dp$covariates, 
                                   dp$dim_covariates, 
                                   lp_type)
  

  if(anyNA(pscore.fitted)){
    stop("Propensity score model have NA values \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }

  # avoid dividing by zero
  PS_MIN = list_control$ps_min
  
  ps.fit <- pmax(pmin(pscore.fitted, 1 - PS_MIN), PS_MIN)
  
  
  #---------------------------------------------------------------------
  # record the call
  call.param <- match.call()
  # record all arguments used in the function
  argu <- mget(names(formals()), sys.frame(sys.nframe()))
  
  # return the output
  ret <- (list(fitted.values = ps.fit,
               bws = bws,
               call.param = call.param,
               argu = argu))
}


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Cross-validated Bandwidth for the OR model
#' Cross-validated bandwidth for the outcome-regression local-polynomial model
#'
#' @param dp Data list with dpost, covariates, and dim_covariates.
#' @param dty Element-wise product of the treatment/time indicators and y.
#' @param dt_index Treatment-group column to fit (-1 for a shared bandwidth).
#' @param n_start Number of random restarts for the bandwidth optimizer.
#' @return Cross-validation bandwidth object returned by cv_or().
locpol.or.bw <- function(dp,
                         dty,
                         dt_index, 
                         n_start
){
  
  #-----------------------------------------------------------------------------
  # initialize the bandwidth parameters 

  if (n_start > 0) {
    rc_bw = runif(n_start-1)
    bwc_init =  c(varc_init, (ubc_init - lbc_init)*rc_bw+lbc_init)
    ru_bw = runif(n_start-1)
    bwu_init =  c(vard_init, (ubd_init - lbd_init)*ru_bw+lbd_init)
    ro_bw = runif(n_start-1)
    bwo_init =  c(vard_init, (ubd_init - lbd_init)*ro_bw+lbd_init)
    
    # put them together
    bw_init = rbind(bwc_init, bwu_init, bwo_init)
    
  } else {
    stop("number of restarts must be a positive integer.")
  }
  
  # estimate the cross-validated bandwidth for nonparametric PS model
  cv.bw <- cv_or(dp$dpost, 
                 dty, 
                 dp$covariates, 
                 dp$dim_covariates,
                 dt_index,
                 n_start,
                 bw_init)
  
  return(cv.bw)
}



#-----------------------------------------------------------------------------
# OR estimates
#' Estimate the four outcome-regression functions by local polynomials
#'
#' Fits the OR for each treatment/time group by least-squares local
#' polynomials, using either a single shared (constrained) bandwidth or
#' separate per-group (unconstrained) bandwidths, selecting by CV when none
#' is supplied.
#'
#' @param dp Data list with dpost, y, covariates, and dim_covariates.
#' @param bws Optional bandwidth(s); NULL selects by CV. A matrix is required
#'   for the unconstrained case.
#' @param list_control List with n_start (restarts) and bw_constrained
#'   (TRUE = shared bandwidth, FALSE = per-group).
#' @return List: fitted.values (n x 4), bws, call.param, argu.
locpol.or.fit <- function(dp,
                          bws = NULL, 
                          list_control = list(n_start = 1,
                                              bw_constrained = FALSE))
                        {
  # calculate dt * y
  dty <- dp$dpost * dp$y
  
 if (list_control$bw_constrained) {
  if(is.null(bws)) {
    bws_list <- locpol.or.bw(dp = dp,
                             dty = dty,
                             dt_index = -1,
                             n_start = list_control$n_start
                             )
    
    bws = bws_list$bw_cv
  } 
  
  # # check if bws has the correct dimension
  # if ( !all(which(bws != 0) == which(dp$dim_covariates != 0))  ) {
  #   stop(paste("must provide at least one bandwidth for each type of covariates") )
  # }
  
  #---------------------------------------------------------------------
  
  # compute the OR estimates by least squares local polynomial regression
  outreg.fitted   <-  locpolfit_or(bws, 
                                   dp$dpost,
                                   dty,
                                   dp$covariates, 
                                   dp$dim_covariates,
                                   -1)
  
  
  if(anyNA(outreg.fitted)){
    stop("Outcome regression model fitted values have NA values
         \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  
  or.fit <- as.matrix(outreg.fitted)
  
} else{
  
  #---------------------------------------------------------------------  
  # check if bandwidth is provided, if not, estimate it by cross validation
  if(is.null(bws)) {
    BW <- matrix(0, nrow = 3, ncol = 4)
  } else {
     if (is.matrix(bws)){
       BW <- bws
     } else {
       stop("bws must bw a matrix object")
     }
  }
   
  
  or.fit <- matrix(0, nrow = length(dp$y), ncol = 4)
  
  for (l in 1:4){
  
 
  if(is.null(bws)) {
    bws_list <- locpol.or.bw(dp   = dp,
                             dty  = dty,
                             dt_index = (l-1),
                             n_start = list_control$n_start
                             )
    
    BW[, l] = bws_list$bw_cv
  } 
  
  # # check if bws has the correct dimension
  #   if ( !all(which(BW[, l] != 0) == which(dp$dim_covariates != 0))  ) {
  #     stop(paste("must provide at least one bandwidth for each type of covariates in group", l) )
  #   }
   
  #---------------------------------------------------------------------
  
  # compute the OR estimates by least squares local polynomial regression
  outreg.fitted   <-  locpolfit_or(BW[, l], 
                                   dp$dpost,
                                   dty,
                                   dp$covariates, 
                                   dp$dim_covariates,
                                   (l-1))
  
  
    if(anyNA(outreg.fitted)){
      stop("Outcome regression model fitted values have NA values
           \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
    }
  
  or.fit[,l] <- as.vector(outreg.fitted)
  
  }
  
  # record optimal bandwidths
  bws <- BW
} 
  
  #---------------------------------------------------------------------
  # record the call
  call.param <- match.call()
  # record all arguments used in the function
  argu <- mget(names(formals()), sys.frame(sys.nframe()))
  
  # return the output
  ret <- (list(fitted.values = or.fit,
               bws = bws,
               call.param = call.param,
               argu = argu))
}



#' tryCatch error handler that reports the message and returns TRUE
#'
#' @param error A condition object.
#' @return TRUE (after printing the error message).
handle_any_error <- function(error) {
  cat("Caught error:", conditionMessage(error), "\n")
  return(TRUE)
}

