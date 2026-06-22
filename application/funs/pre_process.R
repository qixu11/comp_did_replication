#-----------------------------------------------------------------------------
#' Data pre-processing for the DR DiD application.
#'
#' Converts a data frame with named columns into the standardized list format
#' required by the local polynomial estimation functions. Handles:
#'   - Outcome, treatment, and post-period extraction
#'   - Continuous covariate normalization to [0, 1]
#'   - Unordered/ordered categorical covariate binding
#'   - Optional fixed effects
#'   - Validation of group sizes
#'
#' @param yname    Name of outcome variable in data.
#' @param tname    Name of post-period indicator in data.
#' @param dname    Name of treatment indicator in data.
#' @param gname    Name of cluster group variable in data.
#' @param xcnames  Character vector of continuous covariate names.
#' @param xunames  Character vector of unordered categorical covariate names.
#' @param xonames  Character vector of ordered categorical covariate names.
#' @param data     Data frame containing all variables.
#' @param fixed.effect Optional matrix of fixed-effect dummies.
#' @param normalized   Logical; normalize continuous covariates to [0,1].
#' @param weightsname  Optional name of observation weight column.
#' @param boot, boot.type, nboot, inffunc  Metadata stored in output list.
#'
#' @return List with y, d, post, dpost (n x 4), group, covariates (n x k),
#'   dim_covariates (length-3), i.weights, and metadata.
#-----------------------------------------------------------------------------

pre_process_drdid <- function(yname,
                              tname,
                              dname,
                              gname,
                              xcnames = NULL,
                              xunames = NULL,
                              xonames = NULL,
                              data,
                              fixed.effect = NULL, 
                              normalized = TRUE,
                              weightsname = NULL,
                              boot = FALSE,
                              boot.type =  c("weighted", "multiplier"),
                              nboot = NULL,
                              inffunc = FALSE) {


  # set bootstrap type
  boot.type <- boot.type[1]

  # Flag for boot.type
  if (boot){
    if ( (boot.type!="weighted") && (boot.type!="multiplier")) {
      warning("boot.type = ",boot.type,  " is not supported. Using 'weighted'.")
      boot.type <- "weighted"
    }
  }
  
  # Flag for normalized
  if ( (normalized != TRUE) && (normalized != FALSE)) {
    warning("normalized = ",normalized,  " is not supported. Using 'TRUE'.")
    normalized <- TRUE
  }
  
  # make sure dataset is a data.frame
  dta <- data
  # this gets around RStudio's default of reading data as tibble
  if (!all( class(dta) == "data.frame")) {
    #warning("class of data object was not data.frame; converting...")
    dta <- as.data.frame(dta)
  }
  
  
  # Flag for yname
  if ( !is.element(yname, base::colnames(dta))) {
    stop("yname = ",yname,  " could not be found in the data provided.")
  }
  # Flag for tname
  if ( !is.element(tname, base::colnames(dta))) {
    stop("tname = ",tname,  " could not be found in the data provided.")
  }
  # Flag for dname
  if ( !is.element(dname, base::colnames(dta))) {
    stop("dname = ",dname,  " could not be found in the data provided.")
  }
  
  # set weights if null
  base::ifelse(is.null(weightsname), w <- rep(1,nrow(dta)), w <- dta[,weightsname])
  dta$w <- w
  
  
  # make sure time periods are numeric
  if (! (is.numeric(dta[, tname])) ) {
    stop("data[, tname] must be numeric. Please convert it.")
    
  }
  
  #  make sure dname is numeric
  if (! (is.numeric(dta[, dname])) ) {
    stop("data[, dname] must be numeric. Please convert it.")
    
  }
  
  # figure out the time periods
  # list of dates from smallest to largest
  tlist <- unique(dta[,tname])[base::order(unique(dta[,tname]))]
  if ( length(na.omit(tlist))!=2) {
    stop("this package only covers the cases with two time periods (pre and post) and two treatment groups (d=1 if treated at post, and d=0 if not treated in both pre and post).
           See package `did' for the cases with multiple groups and/or multiple time periods.")
  }
  
  # list of groups from smallest to largest
  glist <- unique(dta[,dname])[base::order(unique(dta[,dname]))]
  if ( length(na.omit(glist))!=2) {
    stop("drdid only work for case with two time periods and two treatment groups (d=1 if treated at post, and d=0 if not treated in both pre and post).
         See package `did' for the cases with multiple groups and/or multiple time periods.")
  }
  
  
  
  # Outcome variable will be denoted by y
  dta$y <- dta[, yname]
  # Treatment group will be denoted by d
  dta$d <- dta[, dname]
  # Post dummy will be denoted by post
  dta$post <- dta[, tname]
  # Cluster groups will be denoted by group
  dta$group <- dta[, gname]
  
  # matrix of treatment by time group dummies
  dta$dpost <- cbind(dta$d * dta$post,  
                     dta$d * (1-dta$post), 
                     (1-dta$d) * dta$post,
                     (1-dta$d) * (1-dta$post))
  
  dim_covariates <- rep(0, 3)  # dimensions of the three types of variables
  covariates  <-  matrix(nrow = length(dta$d), ncol = 0)
  # get continuous covariates
  if ( !is.null(xcnames) ){
  cov.cont <- dta[, xcnames]
  dim_covariates[1] <- length(xcnames)
  
  # normalize the continuous variables to lie in [0, 1]
   for (s in 1:dim_covariates[1]){
    cov.cont[, s]  <-  (cov.cont[, s] - min(cov.cont[, s], na.rm = TRUE))/
          (max(cov.cont[, s], na.rm = TRUE) - min(cov.cont[, s], na.rm = TRUE))
   }
  
   covariates <- cbind(covariates, cov.cont)
  }
  
  # get categorical variables 
  if ( !is.null(xunames) ){
    dim_covariates[2] <- length(xunames)
    covariates <- cbind(covariates, dta[, xunames])
  }
 
  if (!is.null(xonames)){
    dim_covariates[3] <- length(xonames)
    covariates <- cbind(covariates,  dta[, xonames])
  }
  
  if ( is.null(xcnames) &is.null(xunames) &is.null(xonames) ){
    stop("include at least one covariate.")
  }
  
  if (!is.null(fixed.effect)){
    covariates <- cbind(covariates,  fixed.effect)
  }
  
  # check against very small groups
  gsize <- stats::aggregate(dta$d, by=list(dta$d), function(x) length(x)/length(tlist))
  
  # how many in each group before give warning
  
  # 5 is just a buffer, could pick something else, but seems to work fine
  reqsize <- sum(dim_covariates) + 5
  
  # which groups to warn about
  gsize <- subset(gsize, gsize$x < reqsize) # x is name of column from gsize
  
  # warn if some groups are small
  if (nrow(gsize) > 0) {
    stop("either treatment or the comparison group in your dataset is very small.
          Inference is not feasible.")
  }
  
  #-----------------------------------------------------------------------------
  # setup data 
    dta.final <- as.data.frame(cbind(y = dta$y, 
                                     d = dta$d, 
                                     post = dta$post, 
                                     dta$dpost, 
                                     group = dta$group,
                                     w = dta$w, 
                                     covariates))
  # remove NAs
    dta.final <- dta.final[stats::complete.cases(dta.final), ]
  
  # return the final data list 
    out <- list(y = dta.final$y,
                d = dta.final$d,
                post = dta.final$post,
                dpost = as.matrix(dta.final[4:7]),
                group = dta.final$group,
                covariates = as.matrix(dta.final[,-c(1:9)]),
                dim_covariates = dim_covariates,
                i.weights = dta.final$w,
                normalized = normalized,
                boot = boot,
                boot.type = boot.type,
                nboot = nboot,
                inffunc = inffunc
                )
  
  
  return(out)
}
