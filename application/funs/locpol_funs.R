#-----------------------------------------------------------------------------
# Local polynomial estimation wrappers for the application.
#
# This file sources the application-specific C++ backend (np_lp.cpp) and
# provides R-level wrappers for GPS and OR bandwidth selection and fitting.
#
# Note: the application uses only continuous (k_c = 2) and unordered (k_u = 9)
# covariates, with no ordered covariates (k_o = 0). The bandwidth vector is
# therefore length 2: (h_continuous, lambda_unordered).
#-----------------------------------------------------------------------------

# Source Rcpp file (application-specific C++ backend)
Rcpp::sourceCpp(file.path(base_dir, "funs/np_lp.cpp"))

# Constants for bandwidth initialization bounds
varc_init <- 0.8
vard_init <- 0.8
lbc_init  <- 0.1
ubc_init  <- 1
lbd_init  <- 0.1
ubd_init  <- 1

#-----------------------------------------------------------------------------
#' Cross-validated bandwidth for the GPS model
#'
#' @param dp Data list from pre_process_drdid().
#' @param cv_method CV loss: "cv.ml" (log-likelihood) or "cv.ls" (least squares).
#' @param lp_type Estimation method: "logit" or "ls".
#' @param n_start Number of Nelder-Mead random restarts.
#' @return List with bw_cv (optimal bandwidths), value, and convergence flag.
locpol.ps.bw <- function(dp,
                         cv_method,
                         lp_type,
                         n_start) {
 
  #mysd <- EssDee(dp$covariates[, 1:dp$dim_covariates[1]])
  
  #---------------------------------------------------------------------  
  # initialize the bandwidth parameters 
  if (n_start > 0){
    
    rc_bw = runif(n_start-1)
    bwc_init =  c(varc_init, (ubc_init - lbc_init)*rc_bw+lbc_init)
    ru_bw = runif(n_start-1)
    bwu_init =  c(vard_init, (ubd_init - lbd_init)*ru_bw+lbd_init)
    # ro_bw = runif(n_start-1)
    # bwo_init =  c(vard_init, (ubd_init - lbd_init)*ro_bw+lbd_init)
    
    # put them together
    bw_init = rbind(bwc_init, bwu_init)
    
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
#' GPS estimation via local polynomial logistic regression
#'
#' Selects bandwidths by LOOCV (if bws = NULL), then fits the local
#' multinomial logit GPS model at each observation.
#'
#' @param dp Data list from pre_process_drdid().
#' @param bws Optional bandwidth vector. If NULL, selected by CV.
#' @param cv_method CV loss for bandwidth selection.
#' @param lp_type Estimation method: "logit" or "ls".
#' @param list_control List with n_start (restarts) and ps_min (trimming bound).
#' @return List with fitted.values (n x 4 GPS matrix) and bws.
locpol.ps.fit <- function(dp,
                          bws = NULL,
                          cv_method = "cv.ml",
                          lp_type = "logit",
                          list_control = list(n_start = 1,
                                              ps_min = 1e-3)) {
  
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
  
  # check if bws has the correct dimension
   if ( !all(which(bws != 0) == which(dp$dim_covariates != 0)) ) {
     stop("must provide at least one bandwidth for each type of covariates in the data.")
   }
  
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
#' Cross-validated bandwidth for the OR model
#'
#' @param dp Data list from pre_process_drdid().
#' @param dty Matrix of dpost * y (n x 4).
#' @param dt_index Group index: 0-3 for unconstrained, -1 for constrained.
#' @param n_start Number of Nelder-Mead random restarts.
#' @return List with bw_cv, value, and convergence flag.
locpol.or.bw <- function(dp,
                         dty,
                         dt_index,
                         n_start) {
  
  #-----------------------------------------------------------------------------
  # initialize the bandwidth parameters 

  if (n_start > 0) {
    rc_bw = runif(n_start-1)
    bwc_init =  c(varc_init, (ubc_init - lbc_init)*rc_bw+lbc_init)
    ru_bw = runif(n_start-1)
    bwu_init =  c(vard_init, (ubd_init - lbd_init)*ru_bw+lbd_init)
    # ro_bw = runif(n_start-1)
    # bwo_init =  c(vard_init, (ubd_init - lbd_init)*ro_bw+lbd_init)
    
    # put them together
    bw_init = rbind(bwc_init, bwu_init)
    
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
#' OR estimation via local polynomial regression
#'
#' Selects bandwidths by LOOCV (if bws = NULL), then fits the local
#' least-squares outcome regression for each (treatment, period) cell.
#'
#' @param dp Data list from pre_process_drdid().
#' @param bws Optional bandwidths. If NULL, selected by CV.
#' @param list_control List with n_start and bw_constrained (shared vs per-group bw).
#' @return List with fitted.values (n x 4 OR matrix) and bws.
locpol.or.fit <- function(dp,
                          bws = NULL,
                          list_control = list(n_start = 2,
                                              bw_constrained = FALSE)) {
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
  
  # check if bws has the correct dimension
  if ( !all(which(bws != 0) == which(dp$dim_covariates != 0))  ) {
    stop(paste("must provide at least one bandwidth for each type of covariates") )
  }
  
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
    BW <- matrix(0, nrow = length(dp$dim_covariates[dp$dim_covariates >0]), ncol = 4)
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
  
  # check if bws has the correct dimension
    if ( !all(which(BW[, l] != 0) == which(dp$dim_covariates != 0))  ) {
      stop(paste("must provide at least one bandwidth for each type of covariates in group", l) )
    }
   
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





