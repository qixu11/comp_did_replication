# install_packages.R
# Install the R packages required by this replication package.
# Run once before reproducing the simulations or the application:  source("install_packages.R")
#
# Note: Rcpp/RcppArmadillo/roptim require a working C++ toolchain
# (macOS: Xcode Command Line Tools; Windows: Rtools; Linux: g++), because the
# local-polynomial backend (np_lp*.cpp) is compiled at run time.

pkgs <- c(
  "Rcpp", "RcppArmadillo", "roptim",   # C++ backend + optimizer
  "foreach", "doParallel", "doRNG",    # parallel Monte Carlo
  "np", "nnet", "MASS",                # density/bandwidth + multinomial logit + ginv
  "glmnet", "haven", "dplyr",          # application: data import / wrangling
  "ggplot2", "reshape2", "patchwork", "scales",  # Figure 1
  "openxlsx",                          # Table 4 (.xlsx)
  "here"                               # project-root-relative file paths
)

to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
if (all(ok)) {
  message("All required packages are installed.")
} else {
  warning("Missing packages: ", paste(pkgs[!ok], collapse = ", "))
}
