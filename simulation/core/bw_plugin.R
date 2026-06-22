#' Plug-in bandwidth for outcome regression (OR) with mixed covariates
#'
#' @param dty Matrix of treatment/time indicators times outcomes.
#' @param dt Matrix of treatment/time indicators.
#' @param x Covariate matrix.
#' @param dim_x Vector of covariate dimensions (continuous, unordered, ordered).
#' @param kernel Kernel name.
#' @param bwmethod Bandwidth selection method for density estimation.
#' @return Numeric bandwidth for OR.
plugin_bw_or <- function(dty, dt, x, dim_x, kernel = "Epanechnikov", bwmethod = "normal-reference") {
  
  # current version works only with continuous covariates (with dim_xc = 2)
  n        <- dim(dt)[1]
  k_c      <- dim_x[1]
  k_u      <- dim_x[2]
  k_o      <- dim_x[3]
  ################################################################
  ## Estimate the density of X_c
  ################################################################  
  x1 <- x[,1]
  x2 <- x[,2]
  x_c <- as.data.frame(cbind(x1,x2))
  
  x_df <- as.data.frame(x)
  for (col in (k_c + 1):(k_c + k_u)) {
    x_df[, col] <- as.factor(x_df[, col])
  }
  
  for (col in (k_c + k_u + 1): sum(dim_x)) {
    x_df[, col] <- as.ordered(x_df[, col])
  }
  
  bwx <- np::npudensbw(dat = x_df, bwmethod = bwmethod)
  fx_fit <- np::npudens(bws = bwx, tdat = x_df, edat = x_df)
  f_x <- fx_fit$dens
  f_x <- pmax(f_x, 1e-5)
  # bwx <- np::npudensbw(dat = x, bwmethod = "cv.ml")
  # fx_fit  <- np::npudens(bws = bwx, tdat = x, edat = x)
  # f_x  <- fx_fit$dens
  
  # Generate all combinations of factors
  combinations <- unique(x[, (k_c + 1) : dim(x)[2]])
  
  # Convert combinations to data frame and name columns
  colnames(combinations) <- colnames(x)[(k_c + 1) : dim(x)[2]]
  
  K <- 9 # n1 + n2 + n3 =  2 + 3 + 4
  
  # Initialize matrices
  gd11_fit <- gd12_fit <- gd22_fit <-  matrix(0, n, 4)
  gddd_fit <- vector("list", n)
  
  for (i_factor in 1:nrow(combinations)) {
    combination <- combinations[i_factor, ]
    
    # Get the indices of rows matching the current combination of factors
    idx_factor <- apply(x[, (k_c+1):dim(x)[2]], 1, function(row) all(row == combination))
    # print(sum(idx_factor))
    # if (i_factor == 57){
    #   ii = 1234
    # }
    # Run polynomial regression for the subset of data
    x_subset <- x[idx_factor, , drop = FALSE]
    dt_subset <- dt[idx_factor, , drop = FALSE]
    dty_subset <- dty[idx_factor, , drop = FALSE]
    
    X <- matrix(design_pol(x_subset[, 1:k_c, drop = FALSE], rep(0, k_c), 3), ncol = (K+1))
    
    or.Pol <- matrix(0, K + 1, 4)
    
    for (j in 1:4) {
      or.Pol[, j] <- MASS::ginv(t(X) %*% (dt_subset[, j] * X)) %*% (t(X) %*% (dty_subset[, j]))
    }
    
    X1dot <- cbind(1, 0, 2 * x_subset[, 1], x_subset[, 2], 0, 3 * x_subset[, 1]^2, 2 * x_subset[, 1] * x_subset[, 2], x_subset[, 2]^2, 0)
    X2dot <- cbind(0, 1, 0, x_subset[, 1], 2 * x_subset[, 2], 0, x_subset[, 1]^2, 2 * x_subset[, 1] * x_subset[, 2], 3 * x_subset[, 2]^2)
    X11dot <- cbind(2, 0, 0, 6 * x_subset[, 1], 2 * x_subset[, 2], 0, 0)
    X12dot <- cbind(0, 1, 0, 0, 2 * x_subset[, 1], 2 * x_subset[, 2], 0)
    X22dot <- cbind(0, 0, 2, 0, 0, 2 * x_subset[, 1], 6 * x_subset[, 2])
    
    # gd11_tmp <- gd12_tmp <- gd22_tmp <- gddd_tmp <- matrix(0, n, 4)
    Const <- matrix(1, dim(x_subset)[1], 4)
    
    
    gd11_tmp <- X11dot %*% or.Pol[4:10, ]
    gd12_tmp <- X12dot %*% or.Pol[4:10, ]
    gd22_tmp <- X22dot %*% or.Pol[4:10, ]
    gddd_tmp <- or.Pol[7:10, ]  # Const %*% 
    
    
    gd11_fit[idx_factor, ] <- gd11_tmp #[idx_factor, ]
    gd12_fit[idx_factor, ] <- gd12_tmp #[idx_factor, ]
    gd22_fit[idx_factor, ] <- gd22_tmp #[idx_factor, ]
    gddd_fit[idx_factor]   <- list(gddd_tmp) #[idx_factor, ]
  }
  
  ################################################################
  ## Compute pilot bandwidths: h.tilde, h.breve
  ################################################################
  
  # First, h.tilde
  Rho = plugin_constants(2, kernel = kernel)
  Rho.b = Rho$rho.b
  Rho.v = Rho$rho.v
  
  ht.num <- rep(0, n)
  for (i in 1:n){
    ht.num[i]  <-  sum((Rho.b %*% gddd_fit[[i]])^2)
  }
  
  ht.const.num   = mean(ht.num)
  
  ht.const.denom = 4 * mean(sqrt(sum( Rho.v^2 ))/f_x)  # multiplied by 4 because we are averaging over the 4 treatment groups
  
  h.tilde = (ht.const.num/ht.const.denom * 2*n / (2*5*6))^(-1/8)
  
  # Next, h.breve
  hb.num <- rep(0, n)
  
  rho = plugin_constants(1, kernel = kernel)
  rho.b = rho$rho.b
  rho.v = rho$rho.v
  
  for (i in 1:n){
    hb.num[i] = sum( (rho.b %*% rbind(gd11_fit[i, ], gd12_fit[i, ], gd22_fit[i, ]))^2 )
  }
  
  hb.const.num   = mean(hb.num)
  hb.const.denom = 4 * rho.v * mean(1/f_x)
  
  h.breve = (hb.const.num/hb.const.denom * 2*2*n / (2*4))^(-1/6)
  
  ################################################################
  ## Estimate m^(q+1)(x) by running local quadratic using h.tilde
  ################################################################
  gdot11_fit  <- gdot12_fit <- gdot22_fit<- matrix(rep(0,4*n), nrow = n)
  K =  5 # k_c + k_c *(k_c+1)/2
  dat_eval = cbind(x_c ,1:n)
  for (i in 1:n){
    X <- as.matrix(design_pol(x_c, x_c[i,], 2))
    
    # wgt.bw <- kernel.weights(x, dat_eval[i,], dim_x, 
    #                          h.tilde, lambda_u = 0, lambda_o = 0,   
    #                          kernel = "Epanechnikov")
    flag = 0
    wgt.bw <-  wgt_kernel_mixed(x, dim_x, c(h.tilde, 0, 0), i-1, flag)
    wgt.bw <- as.numeric(wgt.bw)
    if (flag == 1){
      warning("All weights are zero")
    }
    
    or.fit <- matrix(rep(0,4*(K+1)), nrow = (K+1))
    
    or.fit[,1] = (MASS::ginv(t(X)%*%(wgt.bw *dt[,1]*X)) %*% (t(X) %*% (wgt.bw *dty[,1])))
    or.fit[,2] = (MASS::ginv(t(X)%*%(wgt.bw *dt[,2]*X)) %*% (t(X) %*% (wgt.bw *dty[,2])))
    or.fit[,3] = (MASS::ginv(t(X)%*%(wgt.bw *dt[,3]*X)) %*% (t(X) %*% (wgt.bw *dty[,3])))
    or.fit[,4] = (MASS::ginv(t(X)%*%(wgt.bw *dt[,4]*X)) %*% (t(X) %*% (wgt.bw *dty[,4])))
    
    # compute the derivative estimates (from the lp estimates multiplied by the factorial coefficients)
    for (j in 1:4){
      gdot11_fit[i, j]   = 2*or.fit[4, j]  
      gdot12_fit[i, j]   = 1*or.fit[5, j]   # updated the coefficient
      gdot22_fit[i, j]   = 2*or.fit[6, j] 
    }
  }
  
  h.num <- rep(0, n)
  for (i in 1:n){
    h.num[i] = sum( (rho.b %*% rbind(gdot11_fit[i, ], gdot12_fit[i, ], gdot22_fit[i, ]) )^2 )
  }
  h.const.num   = mean(h.num)
  h.const.denom = 4 * rho.v * mean(1/f_x)
  
  h = (h.const.num/h.const.denom * 2*2*n / (2*4))^(-1/6)
  return(h)
}


#' Plug-in bandwidth for propensity score (PS) with mixed covariates
#'
#' @param dt Matrix of treatment/time indicators.
#' @param x Covariate matrix.
#' @param dim_x Vector of covariate dimensions (continuous, unordered, ordered).
#' @param kernel Kernel name.
#' @param bwmethod Bandwidth selection method for density estimation.
#' @return Numeric bandwidth for PS.
plugin_bw_ps <- function(dt, x, dim_x, kernel = "Epanechnikov", bwmethod = "normal-reference") {
  
  # current version works only with continuous covariates (with dim_xc = 2)
  n        <- dim(dt)[1]
  k_c      <- dim_x[1]
  k_u      <- dim_x[2]
  k_o      <- dim_x[3]
  ################################################################
  ## Estimate the density of X_c
  ################################################################  
  x1 <- x[,1]
  x2 <- x[,2]
  x_c <- as.data.frame(cbind(x1,x2))
  
  x_df <- as.data.frame(x)
  for (col in (k_c + 1):(k_c + k_u)) {
    x_df[, col] <- as.factor(x_df[, col])
  }
  
  for (col in (k_c + k_u + 1): sum(dim_x)) {
    x_df[, col] <- as.ordered(x_df[, col])
  }
  
  bwx <- np::npudensbw(dat = x_df, bwmethod = bwmethod)
  fx_fit <- np::npudens(bws = bwx, tdat = x_df, edat = x_df)
  f_x <- fx_fit$dens
  f_x <- pmax(f_x, 1e-5)
  
  # bwx <- np::npudensbw(dat = x, bwmethod = "cv.ml")
  # fx_fit  <- np::npudens(bws = bwx, tdat = x, edat = x)
  # f_x  <- fx_fit$dens
  
  # Generate all combinations of factors
  combinations <- unique(x[, (k_c + 1) : dim(x)[2]])
  
  # Convert combinations to data frame and name columns
  colnames(combinations) <- colnames(x)[(k_c + 1) : dim(x)[2]]
  
  K <- 9 # n1 + n2 + n3 =  2 + 3 + 4
  
  # Initialize matrices
  gd11_fit <- gd12_fit <- gd22_fit <-  matrix(0, n, 3)
  ps.fit <- matrix(0, n, 4)
  gddd_fit <- vector("list", n)
  
  for (i_factor in 1:nrow(combinations)) {
    combination <- combinations[i_factor, ]
    
    # Get the indices of rows matching the current combination of factors
    idx_factor <- apply(x[, (k_c+1) : dim(x)[2]], 1, function(row) all(row == combination))
    
    # Run polynomial regression for the subset of data
    x_subset <- x[idx_factor, , drop = FALSE]
    dt_subset <- dt[idx_factor, , drop = FALSE]
    
    
    X <- matrix(design_pol(x_subset[, 1:k_c, drop = FALSE], rep(0, k_c), 3), ncol = (K+1))
    
    mnl_fit   <- nnet::multinom(dt_subset ~ X-1, trace = F)
    coeff_mnl <- mnl_fit$wts[-c((1:(K+3)), (2*K+5), (3*K+7))]
    ps.Pol <- matrix(coeff_mnl,nrow=K+1,ncol=3)
    ps.fit.tmp    <- mnl_fit$fitted.values
    
    X1dot <- cbind(1, 0, 2 * x_subset[, 1], x_subset[, 2], 0, 3 * x_subset[, 1]^2, 2 * x_subset[, 1] * x_subset[, 2], x_subset[, 2]^2, 0)
    X2dot <- cbind(0, 1, 0, x_subset[, 1], 2 * x_subset[, 2], 0, x_subset[, 1]^2, 2 * x_subset[, 1] * x_subset[, 2], 3 * x_subset[, 2]^2)
    X11dot <- cbind(2, 0, 0, 6 * x_subset[, 1], 2 * x_subset[, 2], 0, 0)
    X12dot <- cbind(0, 1, 0, 0, 2 * x_subset[, 1], 2 * x_subset[, 2], 0)
    X22dot <- cbind(0, 0, 2, 0, 0, 2 * x_subset[, 1], 6 * x_subset[, 2])
    
    # gd11_tmp <- gd12_tmp <- gd22_tmp <- gddd_tmp <- matrix(0, n, 4)
    #Const <- matrix(1, dim(x_subset)[1], 4)   # 7 to 10
    
    gd11_tmp <- X11dot %*% ps.Pol[4:10, ]
    gd12_tmp <- X12dot %*% ps.Pol[4:10, ]
    gd22_tmp <- X22dot %*% ps.Pol[4:10, ]
    gddd_tmp <- ps.Pol[7:10, ]
    
    # assigning the values to appropriate rows 
    gd11_fit[idx_factor, ] <- gd11_tmp #[idx_factor, ]
    gd12_fit[idx_factor, ] <- gd12_tmp #[idx_factor, ]
    gd22_fit[idx_factor, ] <- gd22_tmp #[idx_factor, ]
    gddd_fit[idx_factor] <- list(gddd_tmp) #[idx_factor, ]
    
    ps.fit[idx_factor,] <-  ps.fit.tmp
    
  }
  
  # trim ps scores
  
  ps.fit = trim_and_normalize(ps.fit)
  
  ################################################################
  ## Compute pilot bandwidths: h.tilde, h.breve
  ################################################################
  
  # First, h.tilde
  Rho = plugin_constants(2, kernel = kernel)
  Rho.b = Rho$rho.b
  Rho.v = Rho$rho.v
  
  ht.denom <- rep(0, n)
  
  
  ht.num <- rep(0, n)
  for (i in 1:n){
    ht.num[i]  <-  sum((Rho.b %*% gddd_fit[[i]])^2)   #  2 x 4  * 4 x 3
  }
  ht.const.num  <- mean(ht.num)  # sum((gddd_fit %*% t(Rho.b))^2)
  
  for (i in 1:n){
    I_x = diag(ps.fit[i, 2:4]) -  outer(ps.fit[i, 2:4], ps.fit[i, 2:4], "*")
    ht.denom[i] = sum(diag(MASS::ginv(I_x) %x% Rho.v))
  }
  
  ht.const.denom = mean(ht.denom/f_x)
  
  h.tilde = (ht.const.num/ht.const.denom * 2 * n / (2*5*6))^(-1/8)
  
  # Next, h.breve
  hb.num <- hb.denom <- rep(0, n)
  rho = plugin_constants(1, kernel = kernel)
  rho.b = rho$rho.b
  rho.v = rho$rho.v
  
  for (i in 1:n){
    I_x = diag(ps.fit[i, 2:4]) -  outer(ps.fit[i, 2:4], ps.fit[i, 2:4], "*")
    hb.num[i] = sum((rho.b %*% rbind(gd11_fit[i, ], gd12_fit[i, ],gd22_fit[i, ]) %*% I_x)^2 )
    hb.denom[i] = sum(diag(I_x))
    
  }
  hb.const.num   = mean(hb.num)
  hb.const.denom = rho.v * mean(hb.denom/f_x)
  
  h.breve = (hb.const.num / hb.const.denom * 2*2*n / (2*4))^(-1/6)
  
  ################################################################
  ## Estimate g^(p+1)(x) using h.tilde
  ################################################################
  gdot11_fit  <- gdot12_fit <- gdot22_fit<- matrix(rep(0,3*n), nrow = n)
  K =  5 # k_c + k_c *(k_c+1)/2
  dat_eval = cbind(x_c ,1:n)
  for (i in 1:n){
    X <- as.matrix(design_pol(x_c, x_c[i,], 2))
    
    
    flag = 0
    wgt.bw <-  wgt_kernel_mixed(x, dim_x, c(h.tilde, 0, 0), i-1, flag)
    wgt.bw <- as.numeric(wgt.bw)
    
    # run global polynomial logit regression
    mnl_fit   <- nnet::multinom(dt ~ X-1, weights = wgt.bw, trace = F)
    coeff_mnl <- mnl_fit$wts[-c((1:(K+3)), (2*K+5), (3*K+7))]
    coeff_mnl <- t(matrix(coeff_mnl,nrow=K+1,ncol=3))
    # compute the propensity score estimates
    for (j in 1:3){
      gdot11_fit[i, j]   = 2*coeff_mnl[j,4]
      gdot12_fit[i, j]   = 1*coeff_mnl[j,5]  # updated the coefficient
      gdot22_fit[i, j]   = 2*coeff_mnl[j,6]
    }
  }
  
  ################################################################
  ## Estimate I(x) using h.breve
  ################################################################
  
  #ps.fit <- t(apply(dat_eval, 1, ps_fit, bw= h.breve, DT = DT, x = x, dim_x = dim_x))
  ks_flag = 0
  bw.breve = c(h.breve, 0, 0)
  ps.fit <- lp_logit_ps_fit(bw.breve, dt, x,  dim_x, ks_flag)
  
  h.num <- h.denom <- rep(0, n)
  for (i in 1:n){
    I_x = diag(ps.fit[i, 2:4]) -  outer(ps.fit[i, 2:4], ps.fit[i, 2:4], "*")
    h.num[i] = sum( (rho.b %*% rbind(gdot11_fit[i, ], gdot12_fit[i, ], gdot22_fit[i, ]) %*% I_x)^2 )
    h.denom[i] = sum(diag(I_x))
  }
  h.const.num   = mean(h.num)
  h.const.denom = rho.v * mean(h.denom/f_x)
  
  h = (h.const.num/h.const.denom * 2*2*n / (2*4))^(-1/6)
  
  return(h)
}

#' Trim and normalize probability vectors
#'
#' @param probs Matrix of probabilities.
#' @param min_prob Lower trimming bound.
#' @param max_prob Upper trimming bound.
#' @return Trimmed and row-normalized probabilities.
trim_and_normalize <- function(probs, min_prob = 0.001, max_prob = 0.999) {
  # Step 1: Trim the probabilities
  probs <- pmin(pmax(probs, min_prob), max_prob)
  
  # Step 2: Normalize the probabilities so they sum up to 1
  row_sums <- rowSums(probs)
  probs <- probs / row_sums
  
  return(probs)
}

#' Plug-in constants for bandwidth calculations
#'
#' @param p Polynomial order.
#' @param kernel Kernel name (currently only Epanechnikov supported).
#' @return List with `rho.b` and `rho.v` constants.
plugin_constants <- function(p, kernel = "Epanechnikov") {
  # for now, this works only when k_c = 2
  k2 = 1/5
  k4 = 3/35
  r0 = 3/5
  r2 = 3/35
  r4 = 1/35
  Q11 = diag(c(1, k2, k2))
  M12 = matrix(c(k2, 0, k2, 
                 0, 0, 0, 
                 0, 0, 0), nrow = 3, byrow = TRUE)
  T11 = r0 * diag(c(1, r2, r2))
  if (p ==1){
    
    rho.b = solve(Q11) %*% M12
    rho.b = rho.b[1,]
    rho.v = solve(Q11) %*% T11 %*% solve(Q11)
    rho.v = rho.v[1,1]
    
  } else if (p == 2){
    Q22 =  matrix(c(k4, 0, k2^2,
                    0, k2^2, 0, 
                    k2^2, 0, k4), nrow = 3, byrow = TRUE)
    Q2  =  rbind(cbind(Q11, M12), cbind(t(M12), Q22))  # 6 x 6 ( (n0 + n1 + n2) = 1 + 2 + 3)
    M23 =  matrix(c(0, 0, 0,0,
                    k4, 0,k2^2, 0,
                    0, k2^2, 0,k4,
                    0, 0, 0, 0, 
                    0, 0, 0, 0,
                    0, 0, 0, 0), nrow = 6, byrow = TRUE)   # 6 x 4 ( (n0 + n1 + n2) x n3) # updated the matrix
    T12 =   matrix(c(r0*r2, 0, r0*r2, 
                     0, 0, 0, 
                     0, 0, 0), nrow = 3, byrow = TRUE)
    T22 =  matrix(c(r0*r4, 0, r2^2,
                    0, r2^2, 0, 
                    r2^2, 0, r0*r4), nrow = 3, byrow = TRUE)
    T2  = rbind(cbind(T11, T12), cbind(t(T12), T22))
    
    rho.b = solve(Q2) %*% M23
    rho.b = rho.b[2:3,]                # 2 x 4   (n1 x n3)
    rho.v = solve(Q2) %*% T2 %*% solve(Q2)
    rho.v = rho.v[2:3,2:3]
  }
  
  return(list(rho.b = rho.b,
              rho.v = rho.v))
}


#' Polynomial design matrix for local polynomial regression
#'
#' @param x Covariate matrix.
#' @param dat_eval Evaluation point.
#' @param deg_lp Polynomial degree.
#' @return Design matrix for local polynomial fit.
design_pol <- function(x, dat_eval, deg_lp) {
  # Ensure x and dat_eval are numeric
  x <- as.data.frame(x)
  dat_eval <- as.numeric(dat_eval)
  
  k <- dim(x)[2]
  x <- as.matrix(x)
  x0 <- matrix(dat_eval, nrow = 1, ncol = k, byrow = TRUE)
  
  # Subtract x0 from each row of x
  X <- sweep(x, 2, x0, "-")
  
  # Generate polynomial terms
  x_1 <- poly(X, 1, raw = TRUE)
  x2 <- poly(X, 2, raw = TRUE)
  x3 <- poly(X, 3, raw = TRUE)
  
  if (deg_lp == 1) {
    design_x <- cbind(1, x_1)
  } else if (deg_lp == 2) {
    x_2 <- x2[, -which(colnames(x2) %in% colnames(x_1)), drop = FALSE]
    design_x <- cbind(1, x_1, x_2)
  } else if (deg_lp == 3) {
    x_2 <- x2[, -which(colnames(x2) %in% colnames(x_1)), drop = FALSE]
    x_3 <- x3[, -which(colnames(x3) %in% colnames(x2)), drop = FALSE]
    design_x <- cbind(1, x_1, x_2, x_3)
  }
  return(design_x)
}




#' Plug-in bandwidth for OR with unconstrained bandwidths
#'
#' @param dty Matrix of treatment/time indicators times outcomes.
#' @param dt Matrix of treatment/time indicators.
#' @param x Covariate matrix.
#' @param dim_x Vector of covariate dimensions (continuous, unordered, ordered).
#' @param kernel Kernel name.
#' @param bwmethod Bandwidth selection method for density estimation.
#' @return Vector of bandwidths for OR.
plugin_bw_or_unconstr <- function(dty, dt, x, dim_x, kernel = "Epanechnikov", bwmethod = "normal-reference") {
  
  # current version works only with continuous covariates (with dim_xc = 2)
  n        <- dim(dt)[1]
  k_c      <- dim_x[1]
  k_u      <- dim_x[2]
  k_o      <- dim_x[3]
  ################################################################
  ## Estimate the density of X_c
  ################################################################  
  x1 <- x[,1]
  x2 <- x[,2]
  x_c <- as.data.frame(cbind(x1,x2))
  
  x_df = as.data.frame(x)
  for (col in (k_c+1) : (k_c+k_u)) {
    x_df[,col] <- as.factor(x_df[,col])
  }
  
  for (col in (k_c+k_u+1) : sum(dim_x)) {
    x_df[,col] <- as.ordered(x_df[,col])
  }
  
  bwx <- np::npudensbw(dat = x_df, bwmethod = bwmethod)
  fx_fit  <- np::npudens(bws = bwx, tdat = x_df, edat = x_df)
  f_x  <- fx_fit$dens
  f_x  <- pmax(f_x, 1e-5)
  
  # Generate all combinations of factors
  combinations <- unique(x[, (k_c+1) : dim(x)[2]])
  
  # Convert combinations to data frame and name columns
  colnames(combinations) <- colnames(x)[(k_c+1) : dim(x)[2]]
  
  K <- 9 # n1 + n2 + n3 =  2 + 3 + 4
  
  # Initialize matrices
  gd11_fit <- gd12_fit <- gd22_fit <-  matrix(0, n, 4)
  gddd_fit <- vector("list", n)
  
  for (i_factor in 1:nrow(combinations)) {
    combination <- combinations[i_factor, ]
    
    # Get the indices of rows matching the current combination of factors
    idx_factor <- apply(x[, (k_c+1):dim(x)[2]], 1, function(row) all(row == combination))
    
    # Run polynomial regression for the subset of data
    x_subset <- x[idx_factor, , drop = FALSE]
    dt_subset <- dt[idx_factor, , drop = FALSE]
    dty_subset <- dty[idx_factor, , drop = FALSE]
    
    
    X <- matrix(design_pol(x_subset[, 1:k_c, drop = FALSE], rep(0, k_c), 3), ncol = (K+1))
    
    or.Pol <- matrix(0, K + 1, 4)
    
    for (j in 1:4) {
      or.Pol[, j] <- MASS::ginv(t(X) %*% (dt_subset[, j] * X)) %*% (t(X) %*% (dty_subset[, j]))
    }
    
    X1dot <- cbind(1, 0, 2 * x_subset[, 1], x_subset[, 2], 0, 3 * x_subset[, 1]^2, 2 * x_subset[, 1] * x_subset[, 2], x_subset[, 2]^2, 0)
    X2dot <- cbind(0, 1, 0, x_subset[, 1], 2 * x_subset[, 2], 0, x_subset[, 1]^2, 2 * x_subset[, 1] * x_subset[, 2], 3 * x_subset[, 2]^2)
    X11dot <- cbind(2, 0, 0, 6 * x_subset[, 1], 2 * x_subset[, 2], 0, 0)
    X12dot <- cbind(0, 1, 0, 0, 2 * x_subset[, 1], 2 * x_subset[, 2], 0)
    X22dot <- cbind(0, 0, 2, 0, 0, 2 * x_subset[, 1], 6 * x_subset[, 2])
    
    # gd11_tmp <- gd12_tmp <- gd22_tmp <- gddd_tmp <- matrix(0, n, 4)
    #Const <- matrix(1, dim(x_subset)[1], 4)
    
    
    gd11_tmp <- X11dot %*% or.Pol[4:10, ]
    gd12_tmp <- X12dot %*% or.Pol[4:10, ]
    gd22_tmp <- X22dot %*% or.Pol[4:10, ]
    gddd_tmp <- or.Pol[7:10, ]
    
    gd11_fit[idx_factor, ] <- gd11_tmp #[idx_factor, ]
    gd12_fit[idx_factor, ] <- gd12_tmp #[idx_factor, ]
    gd22_fit[idx_factor, ] <- gd22_tmp #[idx_factor, ]
    gddd_fit[idx_factor]   <- list(gddd_tmp) #[idx_factor, ]
  }
  
  ################################################################
  ## Compute pilot bandwidths: h.tilde, h.breve
  ################################################################
  
  # First, h.tilde
  Rho = plugin_constants(2, kernel = kernel)
  Rho.b = Rho$rho.b
  Rho.v = Rho$rho.v
  
  ht.num <- matrix(rep(0, 4*n), nrow = n) #rep(0, n)
  for (i in 1:n){
    for (j in 1:4){
      ht.num[i, j]  <-  sum((Rho.b %*% gddd_fit[[i]][, j])^2)
    }
  }
  
  ht.const.num   = colMeans(ht.num)   # 4 x 1
  
  ht.const.denom = mean(sqrt(sum(Rho.v^2))/f_x)  # multiplied by 4 because we are averaging over the 4 treatment groups
  
  h.tilde = (ht.const.num/ht.const.denom * 2*n / (2*5*6))^(-1/8)  # 4 x 1
  
  # Next, h.breve
  hb.num <-  matrix(rep(0, 4*n), nrow = n)  # rep(0, n)
  
  rho = plugin_constants(1, kernel = kernel)
  rho.b = rho$rho.b
  rho.v = rho$rho.v
  
  for (i in 1:n){
    for (j in 1:4){
      hb.num[i, j] = sum( (rho.b %*% rbind(gd11_fit[i, j], gd12_fit[i, j], gd22_fit[i, j]))^2 )
    }
  }
  
  hb.const.num   = colMeans(hb.num)   # 4 x 1
  hb.const.denom = rho.v * mean(1/f_x)
  
  h.breve = (hb.const.num/hb.const.denom * 2*2*n / (2*4))^(-1/6)   # 4 x 1
  
  ################################################################
  ## Estimate m^(q+1)(x) by running local quadratic using h.tilde
  ################################################################
  gdot11_fit  <- gdot12_fit <- gdot22_fit<- matrix(rep(0,4*n), nrow = n)
  K =  5 # k_c + k_c *(k_c+1)/2
  dat_eval = cbind(x_c ,1:n)
  for (i in 1:n){
    X <- as.matrix(design_pol(x_c, x_c[i,], 2))
    
    # wgt.bw <- kernel.weights(x, dat_eval[i,], dim_x, 
    #                          h.tilde, lambda_u = 0, lambda_o = 0,   
    #                          kernel = "Epanechnikov")
    flag = 0
    wgt.bw.1 <-  wgt_kernel_mixed(x, dim_x, c(h.tilde[1], 0, 0), i-1, flag)
    wgt.bw.2 <-  wgt_kernel_mixed(x, dim_x, c(h.tilde[2], 0, 0), i-1, flag)
    wgt.bw.3 <-  wgt_kernel_mixed(x, dim_x, c(h.tilde[3], 0, 0), i-1, flag)
    wgt.bw.4 <-  wgt_kernel_mixed(x, dim_x, c(h.tilde[4], 0, 0), i-1, flag)
    wgt.bw.1 <- as.numeric(wgt.bw.1)
    wgt.bw.2 <- as.numeric(wgt.bw.2)
    wgt.bw.3 <- as.numeric(wgt.bw.3)
    wgt.bw.4 <- as.numeric(wgt.bw.4)
    
    if (flag == 1){
      warning("All weights are zero")
    }
    
    or.fit <- matrix(rep(0,4*(K+1)), nrow = (K+1))
    
    or.fit[,1] = (MASS::ginv(t(X)%*%(wgt.bw.1 *dt[,1]*X)) %*% (t(X) %*% (wgt.bw.1 *dty[,1])))
    or.fit[,2] = (MASS::ginv(t(X)%*%(wgt.bw.2 *dt[,2]*X)) %*% (t(X) %*% (wgt.bw.2 *dty[,2])))
    or.fit[,3] = (MASS::ginv(t(X)%*%(wgt.bw.3 *dt[,3]*X)) %*% (t(X) %*% (wgt.bw.3 *dty[,3])))
    or.fit[,4] = (MASS::ginv(t(X)%*%(wgt.bw.4 *dt[,4]*X)) %*% (t(X) %*% (wgt.bw.4 *dty[,4])))
    
    # compute the derivative estimates (from the lp estimates multiplied by the factorial coefficients)
    gdot11_fit[i, ] = 2*or.fit[4, ]  
    gdot12_fit[i, ] = 1*or.fit[5, ]   # updated the coefficient
    gdot22_fit[i, ] = 2*or.fit[6, ] 
    
  }
  
  h.num <- matrix(rep(0, 4*n), nrow = n)
  
  for (i in 1:n){
    for (j in 1:4){
      h.num[i, j] = sum( (rho.b %*% rbind(gdot11_fit[i, j], gdot12_fit[i, j], gdot22_fit[i, j]) )^2 )
    }
  }
  h.const.num   = colMeans(h.num) # 4 x 1
  h.const.denom = 4 * rho.v * mean(1/f_x)
  
  h = (h.const.num/h.const.denom * 2*2*n / (2*4))^(-1/6)   # 4 x 1
  return(h)
}
