###========================================================================
###========================================================================
### TWFE DiD Estimates -- Sequeira (2016, AER)
###
### Replicates Panel A (TWFE) of Table 4.
###
### Two TWFE specifications are estimated for each of four outcome variables:
###   Spec 1: time-invariant covariates + FE dummies
###   Spec 2: Spec 1 + covariate-by-post interactions (time-varying effects)
###
### Inference: cluster bootstrap (with Mammen weights) at the HS 4-digit level.
###========================================================================
###========================================================================

#-----------------------------------------------------------------------------
# Startup
#-----------------------------------------------------------------------------
rm(list = ls())

# Application root via 'here'. here::i_am() anchors the project root to THIS
# script; base_dir then resolves regardless of where R was launched (run from
# anywhere inside the package folder).
here::i_am("application/funs/main_twfe.R")
base_dir <- here::here("application")

library("dplyr")
library("glmnet")

source(file.path(base_dir, "funs/dr_did.R"))

# Reproducibility
seed1 <- 1234
set.seed(seed1)

#-----------------------------------------------------------------------------
# Data preparation
#-----------------------------------------------------------------------------
data <- haven::read_dta(file.path(base_dir, "data/Bribes_Regression.dta"))

# Outcome variables
y1 <- data$bp
y2 <- data$lba
y3 <- data$lba_value
y3[2783] <- 0       # fix data entry issue
y4 <- data$lba_tonnage

# Treatment and time indicators
d     <- data$tariff_change_2008
post  <- data$post_2008
group <- data$hc_4digits

# Covariate matrix (cols 1-11: base covariates, 12-13: FE bases, 14-15: alt. value vars)
X <- cbind(data$lvalue_tonnage,           #  1: log(value/tonnage)
           data$tariff2007,               #  2: tariff in 2007
           data$differentiated,           #  3
           data$agri,                     #  4
           data$perishable,               #  5
           data$dfs,                      #  6
           data$day_w_arrival,            #  7
           data$monitor,                  #  8
           data$psi,                      #  9
           data$rsa,                      # 10
           data$term,                     # 11
           data$clear_agent,              # 12: clearing agent ID
           data$hc_group,                 # 13: HS chapter group
           log(data$tonnage + 1),         # 14: log(tonnage)
           data$lvalue_shipment_metical)  # 15: log(shipment value in metical)

# Select complete cases and extract working variables
Working_data <- cbind(y1, y2, y3, y4, d, post, group, X)
Working_data <- Working_data[complete.cases(Working_data), ]

y1    <- Working_data[, 1]
y2    <- Working_data[, 2]
y3    <- Working_data[, 3]
y4    <- Working_data[, 4]
y_mat <- cbind(y1, y2, y3, y4)
d     <- Working_data[, 5]
post  <- Working_data[, 6]
group <- Working_data[, 7]

# Fixed-effect dummies for clearing agent and HS chapter group
clear_agent <- Working_data[, 19]
hc_group    <- Working_data[, 20]

ca <- cbind(clear_agent == 2, clear_agent == 3, clear_agent == 4,
            clear_agent == 5, clear_agent == 6, clear_agent == 7,
            clear_agent == 8)

hg <- cbind(hc_group == 2,  hc_group == 3,  hc_group == 4,  hc_group == 5,
            hc_group == 6,  hc_group == 8,  hc_group == 9,  hc_group == 10,
            hc_group == 11, hc_group == 12, hc_group == 13, hc_group == 14,
            hc_group == 15)

#-----------------------------------------------------------------------------
# Inference parameters
#-----------------------------------------------------------------------------
nboot <- 9999
alpha <- c(0.90, 0.95, 0.99)
sequeira_fe <- list()

###========================================================================
### Loop over four outcome variables
###========================================================================

for (i in 1:4) {

  # ---- Spec 1: time-invariant covariates + FE ----
  # Covariates differ by outcome for the value/tonnage variable
  if (i == 1 | i == 2) {
    X_fe <- cbind(Working_data[, 8:18], ca, hg, group)
  } else if (i == 3) {
    X_fe <- cbind(Working_data[, c(9:18, 21)], ca, hg, group)
  } else if (i == 4) {
    X_fe <- cbind(Working_data[, c(9:18, 22)], ca, hg, group)
  }

  twfe1 <- twfe_did_rc(y_mat[, i], post, d, X_fe,
                        boot = FALSE, nboot = NULL, inffunc = TRUE, alpha = alpha)
  att_twfe1    <- twfe1$ATT
  linrep_twfe1 <- twfe1$att.inf.func

  # Cluster bootstrap for Spec 1
  dr.boot   <- boot.did.cluster(linrep_twfe1, group, nboot)
  sse_hat   <- IQR(dr.boot) / (qnorm(0.75) - qnorm(0.25))
  cv        <- quantile(abs(dr.boot / sse_hat), probs = alpha)
  se_twfe1  <- sse_hat
  uci_twfe1 <- att_twfe1 + cv * sse_hat
  lci_twfe1 <- att_twfe1 - cv * sse_hat

  # ---- Spec 2: add covariate-by-post interactions ----
  if (i == 1 | i == 2) {
    X_fe2 <- cbind(X_fe, post * Working_data[, c(8:14, 16)])
  } else if (i == 3) {
    X_fe2 <- cbind(X_fe, post * Working_data[, c(21, 9:14, 16)])
  } else if (i == 4) {
    X_fe2 <- cbind(X_fe, post * Working_data[, c(22, 9:14, 16)])
  }

  twfe2 <- twfe_did_rc(y_mat[, i], post, d, X_fe2,
                        boot = FALSE, nboot = NULL, inffunc = TRUE, alpha = alpha)
  att_twfe2    <- twfe2$ATT
  linrep_twfe2 <- twfe2$att.inf.func

  # Cluster bootstrap for Spec 2
  dr.boot   <- boot.did.cluster(linrep_twfe2, group, nboot)
  sse_hat   <- IQR(dr.boot) / (qnorm(0.75) - qnorm(0.25))
  cv        <- quantile(abs(dr.boot / sse_hat), probs = alpha)
  se_twfe2  <- sse_hat
  uci_twfe2 <- att_twfe2 + cv * sse_hat
  lci_twfe2 <- att_twfe2 - cv * sse_hat

  # ---- Store results ----
  sequeira_fe[[i]] <- list(
    att_twfe1       = att_twfe1,
    se_twfe1        = twfe1$se,
    lci_twfe1       = twfe1$lci,
    uci_twfe1       = twfe1$uci,
    att_twfe2       = att_twfe2,
    se_twfe2        = twfe2$se,
    lci_twfe2       = twfe2$lci,
    uci_twfe2      = twfe2$uci,
    sse_twfe1     = se_twfe1,
    slci_twfe1    = lci_twfe1,
    suci_twfe1    = uci_twfe1,
    sse_twfe2     = se_twfe2,
    slci_twfe2    = lci_twfe2,
    suci_twfe2   = uci_twfe2
  )
  
}

# Save results
save.image(file.path(base_dir, "results/result_sequeira_twfe.RData"))
