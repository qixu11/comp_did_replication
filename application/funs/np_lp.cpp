/* This file contains functions for local polynomial estimation */
#include <cmath>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]

using namespace roptim;
using namespace arma;

const int     MAXITER   = 1000;
const double  PS_MIN    = 1e-05;       
const double  FTOL1     = 1e-02; 
const double  FTOL2     = 1e-07; 
const double  DBL_SMALL = 1e-05;  
const double  DBL_LARGE = 1e8;
const double  DBL_UPPER = DBL_MAX;

////==============================================================================================
//  AUXILLIARY FUNCTIONS
////==============================================================================================
// [[Rcpp::export]]
vec wgt_kernel_mixed( mat  &x,
                      vec  &dim_x,
                      vec  &bw,
                      int  &idx,
                      bool &flag) {
  
  size_t num_obs    = x.n_rows;
  mat  x_c = x.cols(0, dim_x(0)-1);
  mat  x_u = x.cols(dim_x(0), dim_x(0)+dim_x(1)-1); 
  // mat  x_o = x.cols(dim_x(0)+dim_x(1), dim_x(0)+dim_x(1)+dim_x(2)-1);
  
  int    supp_h = 1;
  vec wgt =  ones< vec>(num_obs);
  
  for (int i = 0; i < num_obs; ++i){
    for (int l = 0; l < dim_x(0); ++l){
      supp_h = ( (x_c(i, l) <= (x_c(idx, l) + bw(0))) * (x_c(i, l) >= (x_c(idx,l) - bw(0))) );
      wgt(i) = wgt(i) * 0.75 * (1- std::pow((x_c(i, l) - x_c(idx, l))/bw(0), 2)) / bw(0) * supp_h;
    }
    
    for (int m = 0; m < dim_x(1); ++m){
      wgt(i) = wgt(i) * (1 + (bw(1) - 1) * std::fabs(x_u(i, m)-x_u(idx, m) ) );
    }
    

  } 
  
  // leave idx'th observation out
  wgt(idx) = 0.0;
  
  if (accu(wgt != 0.0) < 10){  
    flag = 1;
  }
  return wgt;
}


////==============================================================================================
vec mnlmin(size_t    n, 
           vec &b0, 
           vec &_wgt_bw,
           mat &_dpost,
           mat &_X
)
{
  int    maxit  = 10000;
  double abstol = -INFINITY;
  double reltol = sqrt(2.220446e-16);
  double stepredn	= 0.2;
  double acctol	  = 0.0001;
  double reltest  = 10.0;
  
  double Fmin;
  bool   accpoint, enough, fail;
  int    count, funcount, gradcount;
  double f, gradproj;
  int    i, j, ilast, iter = 0;
  double s, steplength;
  double D1, D2;
  
  size_t _num_obs = _X.n_rows;
  size_t _npar = n/3;
  mat _wgtX    = _wgt_bw % _X.each_col();
  vec _gamma   =  zeros< vec>(3*_npar);
  mat _xg      =  zeros< mat>(_num_obs, 3);
  mat _exg     =  ones< mat>(_num_obs, 3);
  vec _sexg    =  ones< mat>(_num_obs) * 3;
  
  vec b = b0;
  vec g(n);
  vec t(n);
  vec X(n);
  vec c(n);
  mat B(n,n);
  
  /// =============================================================  
  // initial evaluation of f
  // update common components
  _xg     = _X *  reshape(b, _npar, 3);   // nobs x 3
  _exg    =  exp(_xg);
  _sexg   =  sum(_exg, 1);
  // compute value of the objective function 
  f = -  sum(( sum(_dpost.cols(1,3) % _xg, 1) -  log( 1 + _sexg) ) % _wgt_bw);
  if ( isnan(f) | isinf(f)  ) f =  DBL_LARGE;
  /// =============================================================
  Fmin = f;
  funcount = gradcount = 1;
  // initial update of the gradient:   fmingr(n0, b, g, ex);  
  /// =============================================================
  mat dev = _dpost.cols(1,3) - _exg.each_col()/(1 + _sexg);  // nobs x 3
  g.subvec(0, _npar-1)         =  - trans( sum(dev.col(0) % _wgtX.each_col(), 0));
  g.subvec(_npar, 2*_npar-1)   =  - trans( sum(dev.col(1) % _wgtX.each_col(), 0));
  g.subvec(2*_npar, 3*_npar-1) =  - trans( sum(dev.col(2) % _wgtX.each_col(), 0));
  /// =============================================================
  iter++;
  ilast = gradcount;
  
  do {
    if (ilast == gradcount) {
      B  =  eye(n, n);
    }
    X = b;
    c = g;
    gradproj = 0.0;
    
    for (i = 0; i < n; i++) {
      s = 0.0;
      for (j = 0; j <= i; j++) s -= B(i, j) * g(j);
      for (j = i + 1; j < n; j++) s -= B(j, i) * g(j);
      t(i) = s;
      gradproj += s * g(i);
    }
    
    if (gradproj < 0.0) {	/* search direction is downhill */
steplength = 1.0;
      accpoint = FALSE;
      do {
        count = 0;
        for (i = 0; i < n; i++) {
          b(i) = X(i) + steplength * t(i);
          if (reltest + X(i) == reltest + b(i)) /* no change */
count++;
        }
        if (count < n) {
          // update objective function
          /// =============================================================
          // update common components
          _xg     = _X *  reshape(b, _npar, 3);   // nobs x 3
          _exg    =  exp(_xg);
          _sexg   =  sum(_exg, 1);
          // compute value of the objective function 
          f = -  sum(( sum(_dpost.cols(1,3) % _xg, 1) -  log( 1 + _sexg) ) % _wgt_bw);
          if ( isnan(f) ) f =  DBL_LARGE;
          /// =============================================================
          funcount++;
          accpoint = std::isfinite(f) &&
            (f <= Fmin + gradproj * steplength * acctol);
          if (!accpoint) {
            steplength *= stepredn;
          }
        }
      } while (!(count == n || accpoint));
      enough = (f > abstol) &&
        std::fabs(f - Fmin) > reltol * (std::fabs(Fmin) + reltol);
      /* stop if value if small or if relative change is low */
      if (!enough) {
        count = n;
        Fmin = f;
      }
      if (count < n) {/* making progress */
      Fmin = f;
        //  update the gradient
        /// =============================================================
        dev = _dpost.cols(1,3) - _exg.each_col()/(1 + _sexg);  // nobs x 3
        g.subvec(0, _npar-1)         =  - trans( sum(dev.col(0) % _wgtX.each_col(), 0));
        g.subvec(_npar, 2*_npar-1)   =  - trans( sum(dev.col(1) % _wgtX.each_col(), 0));
        g.subvec(2*_npar, 3*_npar-1) =  - trans( sum(dev.col(2) % _wgtX.each_col(), 0));
        /// =============================================================
        gradcount++;
        iter++;
        D1 = 0.0;
        for (i = 0; i < n; i++) {
          t(i) = steplength * t(i);
          c(i) = g(i) - c(i);
          D1 += t(i) * c(i);
        }
        if (D1 > 0) {
          D2 = 0.0;
          for (i = 0; i < n; i++) {
            s = 0.0;
            for (j = 0; j <= i; j++)
              s += B(i,j) * c(j);
            for (j = i + 1; j < n; j++)
              s += B(j,i) * c(j);
            X(i)= s;
            D2 += s * c(i);
          }
          D2 = 1.0 + D2 / D1;
          for (i = 0; i < n; i++) {
            for (j = 0; j <= i; j++)
              B(i,j) += (D2 * t(i) * t(j)
                           - X(i) * t(j) - t(i) * X(j)) / D1;
          }
        } else {	/* D1 < 0 */
      ilast = gradcount;
        }
      } else {	/* no progress */
      if (ilast < gradcount) {
        count = 0;
        ilast = gradcount;
      }
      }
    } else {		/* uphill search */
      count = 0;
      if (ilast == gradcount) count = n;
      else ilast = gradcount;
      /* Resets unless has just been reset */
    }
    if (iter >= maxit) break;
    if (gradcount - ilast > 2 * n)
      ilast = gradcount;	/* periodic restart */
  } while (count != n || ilast != gradcount);
  
  fail = (iter >= maxit);
  
  return b;
  
}


////==============================================================================================
//  GPS FUNCTIONS
////==============================================================================================
// [[Rcpp::export]]
mat lp_logit_ps_fit( vec  &bw, 
                     mat  &dpost,  
                     mat  &cov, 
                     vec  &dim_x,
                     bool &ks_flag){
  size_t num_obs = cov.n_rows;
  size_t npar = dim_x(0)+1;
  size_t Npar = 3*npar;
  
  mat  X(num_obs, npar,  fill::ones);
  mat  ps_fit(num_obs, 4,  fill::zeros);
  vec  gamma_0(Npar,  fill::zeros);
  vec  flag_vec(num_obs,  fill::zeros);
  X.cols(1, dim_x(0)) = cov.cols(0, dim_x(0)-1);
// omp_set_num_threads(8);
// #pragma omp parallel for shared(cov, X, dim_x, bw, dpost, ps_fit, flag_vec, npar, num_obs, gamma_0) schedule(dynamic) 
  for (int j = 0; j < num_obs; ++j){
    // get kernel weights
    bool flag = 0;
    vec wgt_bw    =  wgt_kernel_mixed(cov, dim_x, bw, j, flag);
    flag_vec(j)   =  flag;
    // obtain multinomial logit estimates
    vec gamma_opt =  mnlmin(Npar, gamma_0, wgt_bw, dpost, X);
    mat Gamma     =  reshape(gamma_opt, npar, 3);
    vec egx_fit   =  trans( exp(Gamma.row(0)+cov.submat(j,0,j,dim_x(0)-1)*Gamma.rows(1,dim_x(0)))); 
    double    egx_sum   =   sum(egx_fit);
    // assign predicted values
    ps_fit(j, 1) = egx_fit(0)/(1 + egx_sum);
    ps_fit(j, 2) = egx_fit(1)/(1 + egx_sum);
    ps_fit(j, 3) = egx_fit(2)/(1 + egx_sum);
    ps_fit(j, 0) = 1/(1+egx_sum);
  }
  
  ks_flag = any(flag_vec == 1); 
  
  return ps_fit;
}

////===============================================================================================
// [[Rcpp::export]]
mat lp_ls_ps_fit( vec &bw,
                  mat &dpost,
                  mat &cov,
                  vec &dim_x,
                  bool &ks_flag){
  size_t     num_obs = cov.n_rows;
  mat  X(num_obs, dim_x(0)+1,  fill::ones);
  mat  cov_diff(num_obs, cov.n_cols);
  mat  ps_fit(num_obs, 4);
  vec  treat_post(num_obs);
  vec  treat_pre(num_obs);
  vec  cont_post(num_obs);
  vec  cont_pre(num_obs);
  vec  flag_vec(num_obs,  fill::zeros);
  for (int j = 0; j < num_obs; ++j){
    bool flag = 0;
    vec  wgt_bw = wgt_kernel_mixed(cov, dim_x, bw, j, flag);
    flag_vec(j)   =  flag;
    cov_diff = cov.each_row() - cov.row(j);
    X.cols(1, dim_x(0)) = cov_diff.cols(0, dim_x(0)-1);
    
    treat_post =  pinv(X.t() * (wgt_bw %  X.each_col()), DBL_SMALL, "std") *  (X.t() *  (wgt_bw % dpost.col(0)));
    treat_pre  =  pinv(X.t() * (wgt_bw %  X.each_col()), DBL_SMALL, "std") *  (X.t() *  (wgt_bw % dpost.col(1)));
    cont_post  =  pinv(X.t() * (wgt_bw %  X.each_col()), DBL_SMALL, "std") *  (X.t() *  (wgt_bw % dpost.col(2)));
    cont_pre   =  pinv(X.t() * (wgt_bw %  X.each_col()), DBL_SMALL, "std") *  (X.t() *  (wgt_bw % dpost.col(3)));
    
    ps_fit(j, 0) = treat_post(0);
    ps_fit(j, 1) = treat_pre(0);
    ps_fit(j, 2) = cont_post(0);
    ps_fit(j, 3) = cont_pre(0);
  }
  
  ks_flag = any(flag_vec == 1); 
  
  return ps_fit;
}

////===============================================================================================
class LocPol_ps : public Functor {
public:
  LocPol_ps(const  mat & dpost,
            const  mat & cov,
            const  vec & dim_x,
            const std::string &bw_method,
            const std::string &lp_method) :
  _dpost(dpost),  _cov(cov), _dim_x(dim_x), _bw_method(bw_method), _lp_method(lp_method){}
  double operator()(const  vec &bw) override {
    
    size_t    num_obs = _cov.n_rows;
    mat ps_fit(num_obs, 4);
    mat ps_err(num_obs, 4);
    vec BW = bw;
    double    obj;
    int   num_c = 1;
    // check if bandwidth lies in the unit interval
    bool  CHECK_BW_MAX = (min(BW.subvec(0, num_c-1)) > DBL_UPPER) || ((BW.n_elem > num_c)?  (max(BW.subvec(num_c, BW.n_elem - 1))) > 1 : 0 ); 
    bool  CHECK_BW_MIN = min(BW) < 0.0;
    bool  ks_flag = 0;
    
    if (CHECK_BW_MAX || CHECK_BW_MIN){
      obj = DBL_LARGE;
    } else{
      // if boundary is not violated, proceed to compute mse
      if (_lp_method == "logit"){
        ps_fit = lp_logit_ps_fit(BW,  _dpost, _cov, _dim_x, ks_flag);
      }else{
        ps_fit = lp_ls_ps_fit(BW, _dpost, _cov, _dim_x, ks_flag);
        ps_fit.elem(  find(ps_fit <= PS_MIN) ).fill(PS_MIN);
      }
      
      
      if (_bw_method == "cv.ml"){
        ps_fit.elem(  find(ps_fit <= PS_MIN) ).fill(PS_MIN);
        ps_err = _dpost %  log(ps_fit);
        obj = (ks_flag) ?  DBL_LARGE : - mean( mean(ps_err));
        
      }else{
        ps_err = _dpost - ps_fit;
        obj = (ks_flag) ?  DBL_LARGE : accu( square(ps_err))/num_obs; 
      }
      
      if ( isnan(obj) | isinf(obj) ){
        obj = DBL_LARGE;
      }
      
    }
    return obj;
  }
  
private:
  mat _dpost, _cov;
  vec _dim_x;
  std::string _bw_method, _lp_method;
};

////===============================================================================================
// [[Rcpp::export]]
Rcpp:: List cv_ps( mat         &dpost,
                   mat         &cov,
                   vec         &dim_x,
                   std::string &bw_method,
                   std::string &lp_method,
                   size_t      &n_start,
                   mat         &bw_init){
  
  //size_t    num_obs = cov.n_rows;
  
  // bandwith and objective values
  vec       bw, bw_best;
  double    fcv_best = 0.0;
  
  LocPol_ps OBJ_cv_ps(dpost, cov, dim_x, bw_method, lp_method);
  Roptim<LocPol_ps> opt("Nelder-Mead");
  
  // execute Nelder-Mead optimization using various starting points
  for (int iter_start = 0; iter_start < n_start; ++iter_start){
    // initialize the bandwidth parameter
    bw = bw_init.col(iter_start);
    
  
    opt.control.maxit  = MAXITER;
    opt.control.reltol = FTOL1;
    opt.minimize(OBJ_cv_ps, bw);

    // record the temporary minimizer
    if (iter_start == 0 || opt.value() < fcv_best ){
      bw_best = bw; 
      fcv_best = opt.value();
    }
  }
  
  
  // re-run minimization with smaller error tolerance
  opt.control.reltol = FTOL2;
  opt.minimize(OBJ_cv_ps, bw_best);

  
  // return the results
  return Rcpp::List::create(Rcpp::Named("bw_cv") = bw_best,
                            Rcpp::Named("value") = opt.value(),
                            Rcpp::Named("convergence") = opt.convergence());
}


////===============================================================================================
// [[Rcpp::export]]
mat locpolfit_ps( vec   &bws, 
                  mat   &dpost, 
                  mat   &covariates, 
                  vec   &dim_x, 
                  std::string &lp_method){
  
  size_t     num_obs = covariates.n_rows;
  mat  ps_fit(num_obs, 4);
  bool ks_flag= 0; 
  if (lp_method == "logit"){
    ps_fit = lp_logit_ps_fit(bws, dpost, covariates, dim_x, ks_flag);
  }else if (lp_method == "ls"){
    ps_fit = lp_ls_ps_fit(bws, dpost, covariates, dim_x, ks_flag);
  }
  
  return ps_fit;
}


////==============================================================================================
//  OR FUNCTIONS
////==============================================================================================
// [[Rcpp::export]]
mat lp_constr_or_fit(  vec &bw,
                       mat &dpost,
                       mat &dty,
                       mat &cov,
                       vec &dim_x,
                       bool &ks_flag){
  
  size_t     num_obs = cov.n_rows;
  mat  X(num_obs, dim_x(0)+1, fill::ones);
  mat  cov_diff(num_obs, cov.n_cols);
  mat  or_fit(num_obs, 4);
  vec  treat_post(num_obs), treat_pre(num_obs), cont_post(num_obs), cont_pre(num_obs);
  vec  y11_wgt(num_obs), y10_wgt(num_obs), y01_wgt(num_obs), y00_wgt(num_obs);
  vec  wgt_bw(num_obs);
  vec  flag_vec(num_obs, fill::zeros);
  
  for (int j = 0; j < num_obs; ++j){
    bool flag = 0;
    wgt_bw = wgt_kernel_mixed(cov, dim_x, bw, j, flag);
    
    cov_diff = cov.each_row() - cov.row(j);
    X.cols(1, dim_x(0)) = cov_diff.cols(0, dim_x(0)-1);
    
    y11_wgt = wgt_bw % dty.col(0);
    y10_wgt = wgt_bw % dty.col(1);
    y01_wgt = wgt_bw % dty.col(2);
    y00_wgt = wgt_bw % dty.col(3);
    
    treat_post = pinv(X.t() * (wgt_bw % dpost.col(0) % X.each_col()), DBL_SMALL, "std") *  (X.t() *  y11_wgt);
    treat_pre  = pinv(X.t() * (wgt_bw % dpost.col(1) % X.each_col()), DBL_SMALL, "std") *  (X.t() *  y10_wgt);
    cont_post  = pinv(X.t() * (wgt_bw % dpost.col(2) % X.each_col()), DBL_SMALL, "std") *  (X.t() *  y01_wgt);
    cont_pre   = pinv(X.t() * (wgt_bw % dpost.col(3) % X.each_col()), DBL_SMALL, "std") *  (X.t() *  y00_wgt);
    
    or_fit(j, 0) =  treat_post(0);
    or_fit(j, 1) =  treat_pre(0);
    or_fit(j, 2) =  cont_post(0);
    or_fit(j, 3) =  cont_pre(0);
    
    // set flag for group specific kernel smoothing
    flag  = all(arma::abs(y11_wgt) == 0) || all(arma::abs(y10_wgt) == 0) || all(arma::abs(y01_wgt) == 0) || all(arma::abs(y00_wgt) == 0);
    flag_vec(j) = flag;
    
  }
  ks_flag = any(flag_vec == 1);
  return or_fit;
}



////===============================================================================================
// [[Rcpp::export]]
vec lp_unconstr_or_fit( vec& bw, 
                        vec& dpost, 
                        vec& dty, 
                        mat& cov, 
                        vec& dim_x,
                        bool& ks_flag){
  
  size_t     num_obs = cov.n_rows;
  mat  X(num_obs, dim_x(0)+1,  fill::ones);
  mat  cov_diff(num_obs, cov.n_cols);
  vec  ls_fit;
  vec  or_fit(num_obs);
  vec  y_wgt, wgt_bw;
  vec  flag_vec(num_obs, fill::zeros);
  
  for (int j = 0; j < num_obs; ++j){
    bool flag = 0;
    wgt_bw = wgt_kernel_mixed(cov, dim_x, bw, j, flag);
    cov_diff = cov.each_row() - cov.row(j);
    X.cols(1, dim_x(0)) = cov_diff.cols(0, dim_x(0)-1);
    y_wgt = wgt_bw % dty;
    
    // set flag for group specific kernel smoothing
    flag  = all(arma::abs(y_wgt) == 0);
    ls_fit =  pinv(X.t() * (wgt_bw % dpost % X.each_col()), DBL_SMALL, "std") *  (X.t() * y_wgt);
    or_fit(j) = flag ? DBL_LARGE : ls_fit(0);
    flag_vec(j) = flag;
  }
  
  ks_flag = any(flag_vec == 1);
  return or_fit;
}





////===============================================================================================
class LocPol_or : public Functor {
public:
  LocPol_or(const  mat &dpost, 
            const  mat &dty, 
            const  mat &cov, 
            const  vec &dim_x,
            const  int &or_index) :
  _or_index(or_index), _dpost(dpost), _dty(dty), _cov(cov),  _dim_x(dim_x){}
  double operator()(const  vec &bw) override {
    
    size_t    num_obs = _cov.n_rows;
    size_t    num_c = 1;
    vec       BW = bw;
    double    obj = DBL_LARGE;
    
    // check if bandwidths are admissible
    bool CHECK_BW_MAX = (min(BW.subvec(0, num_c-1)) > DBL_UPPER) || ((BW.n_elem > num_c)?  (max(BW.subvec(num_c, BW.n_elem - 1))) > 1 : 0 ); 
    bool CHECK_BW_MIN = (min(BW) < 0.0);
    bool ks_flag = 0;
    
    if (CHECK_BW_MAX || CHECK_BW_MIN){
      obj = DBL_LARGE;
    } else {
      if ( _or_index== 0 || _or_index == 1 || _or_index == 2 || _or_index == 3 ){
        // calculate mse for a specific treatment group
        vec dt_vec  = _dpost.col(_or_index);
        vec dty_vec = _dty.col(_or_index);
        
        vec or_fit =  lp_unconstr_or_fit(BW, 
                                         dt_vec, 
                                         dty_vec, 
                                         _cov, 
                                         _dim_x,
                                         ks_flag);
        vec or_err = dty_vec - dt_vec % or_fit; 
        
        obj =  (ks_flag) ? DBL_LARGE : dot(or_err, or_err)/num_obs;
        
      } else {
        // calculate mse for all four groups
        mat or_fit =  lp_constr_or_fit(BW, 
                                       _dpost, 
                                       _dty,
                                       _cov, 
                                       _dim_x,
                                       ks_flag);
        mat or_err = _dty - _dpost % or_fit;
        
        obj =  (ks_flag) ? DBL_LARGE : accu(square(or_err))/num_obs; 
        
      } 
      
      // check if obj is finite 
      if ( isnan(obj) || isinf(obj) ){
        obj =  DBL_LARGE;
      }
    }
    return obj;
  }
private:
  int _or_index;
  mat _dpost, _dty, _cov;
  vec _dim_x;
};

////===============================================================================================
// [[Rcpp::export]]
Rcpp:: List cv_or( mat      &dpost, 
                   mat      &dty, 
                   mat      &cov, 
                   vec      &dim_x,
                   int      &or_index,
                   size_t   &n_start,
                   mat      &bw_init
){
  
  //size_t    num_obs = cov.n_rows;
  
  // bandwith and objective values
  vec       bw, bw_best;
  double    fcv_best = 0.0;

  
  LocPol_or OBJ_cv_or(dpost, dty, cov, dim_x, or_index);
  Roptim<LocPol_or> opt("Nelder-Mead");
  
  // execute Nelder-Mead optimization using various starting points
  for (int iter_start = 0; iter_start < n_start; ++iter_start){
    // initialize the bandwidth parameter
    bw = bw_init.col(iter_start);
    
    // run roptim
    opt.control.maxit  = MAXITER;
    opt.control.reltol = FTOL1;
    opt.minimize(OBJ_cv_or, bw);
    
    // record the temporary minimizer
    if (iter_start == 0 || opt.value() < fcv_best ){
      bw_best = bw; 
      fcv_best = opt.value();

    }
  }
  
  // re-run minimization with smaller error tolerance
  opt.control.reltol = FTOL2;
  opt.minimize(OBJ_cv_or, bw_best);


  // return the final cv bandwidth 
  return Rcpp::List::create(Rcpp::Named("bw_cv") = bw_best,
                            Rcpp::Named("value") = opt.value(),
                            Rcpp::Named("convergence") = opt.convergence());
  
}


////===============================================================================================
// [[Rcpp::export]]
mat locpolfit_or( vec &bws,
                  mat &dpost, 
                  mat &dty, 
                  mat &covariates, 
                  vec &dim_x,
                  int &or_index){
  
  size_t num_obs = covariates.n_rows;
  bool  ks_flag =0;
  if (or_index == 0 || or_index == 1 || or_index == 2 || or_index == 3){
    mat or_fit(num_obs, 1);
    vec dt_vec = dpost.col(or_index);
    vec dty_vec = dty.col(or_index);
    or_fit.col(0)   = lp_unconstr_or_fit(bws, 
               dt_vec, 
               dty_vec, 
               covariates, 
               dim_x,
               ks_flag);
    return or_fit;
  } else {
    mat or_fit(num_obs, 4);
    or_fit   = lp_constr_or_fit(bws, 
                                dpost, 
                                dty, 
                                covariates, 
                                dim_x,
                                ks_flag);    
    return or_fit;
  }
  
}


