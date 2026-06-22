#-----------------------------------------------------------------------------
# This file gathers results from various job batches
# and creates per-DGP summary tables matching the paper layout (Tables 1 & 2).
#
# Each job-level RData file contains a matrix `mc` with 500 rows (MC reps)
# and 60 columns. We stack all 10 jobs (5000 rows total) per DGP, compute
# summary statistics, and output one CSV per DGP with three panels:
#   Panel A: TWFE estimators
#   Panel B: Nonparametric DR DiD estimators
#   Panel C: Hausman-type test
#-----------------------------------------------------------------------------
rm(list = ls())

#-----------------------------------------------------------------------------
# Parameters
#-----------------------------------------------------------------------------
n_job <- 10        # number of job batches
n     <- 1000      # sample size
nrep  <- 500       # MC replications per job
dgps  <- c(1, 2)   # DGP indices

# DGP descriptions (for table headers)
dgp_labels <- c(
  "1" = "Non-Stationary Covariate Distribution",
  "2" = "Stationary Covariate Distribution"
)

# Paths
address0 <- here::here("simulation", "Main Simulation Results (Table1&2)", "results")  # via 'here': from the project root
out_dir  <- here::here("simulation", "Main Simulation Results (Table1&2)", "tables")   # via 'here': from the project root

#-----------------------------------------------------------------------------
# Column-index lookup: (estimator_col, asy_var_col, coverage_col, CI_length_col)
#   DR estimators
#     nstnr_ml: cols 3, 4, 5, 6
#     nstnr_ls: cols 7, 8, 9, 10
#     stnr_ml:  cols 11, 12, 13, 14
#     stnr_ls:  cols 15, 16, 17, 18
#   TWFE estimators
#     twfe_base:      cols 52, 53, 54, 55
#     twfe_saturated: cols 56, 57, 58, 59
#   Semiparametric efficiency bound: col 35
#   True ATT: col 1
#-----------------------------------------------------------------------------

# Panel A: TWFE estimators
twfe_idx <- list(
  twfe_base      = c(52, 53, 54, 55),
  twfe_saturated = c(56, 57, 58, 59)
)

# Panel B: Nonparametric DR DiD estimators
dr_idx <- list(
  nstnr_ml = c(3, 4, 5, 6),
  nstnr_ls = c(7, 8, 9, 10),
  stnr_ml  = c(11, 12, 13, 14),
  stnr_ls  = c(15, 16, 17, 18)
)

# Panel C: Hausman-type test (test_stat, rej_10, rej_05, rej_01)
test_idx <- list(
  stb_ml = c(36, 40, 41, 42),
  stb_ls = c(37, 43, 44, 45)
)

#-----------------------------------------------------------------------------
# Helper: build estimation rows from a named list of column indices
#-----------------------------------------------------------------------------
build_est_rows <- function(idx_list, mean.mc, median.bias.mc, rmse.mc) {
  rows <- NULL
  for (m in names(idx_list)) {
    idx <- idx_list[[m]]
    row <- data.frame(
      Method    = m,
      Avg.Bias  = round(mean.mc[idx[1]] - mean.mc[1], 3),
      Med.Bias  = round(median.bias.mc[idx[1]], 3),
      RMSE      = round(rmse.mc[idx[1]], 3),
      Asy.Var   = round(mean.mc[idx[2]], 3),
      Cover.    = round(mean.mc[idx[3]], 3),
      CIL       = round(mean.mc[idx[4]], 3),
      stringsAsFactors = FALSE
    )
    rows <- rbind(rows, row)
  }
  rows
}

#-----------------------------------------------------------------------------
# Helper: build test rows from a named list of column indices
#-----------------------------------------------------------------------------
build_test_rows <- function(idx_list, mean.mc) {
  rows <- NULL
  for (m in names(idx_list)) {
    idx <- idx_list[[m]]
    row <- data.frame(
      Method                   = m,
      Avg.Test.Stats.          = round(mean.mc[idx[1]], 3),
      Emp.Rej.Freq.0.10        = round(mean.mc[idx[2]], 3),
      Emp.Rej.Freq.0.05        = round(mean.mc[idx[3]], 3),
      Emp.Rej.Freq.0.01        = round(mean.mc[idx[4]], 3),
      stringsAsFactors = FALSE
    )
    rows <- rbind(rows, row)
  }
  rows
}

#-----------------------------------------------------------------------------
# Loop over DGPs — one table per DGP
#-----------------------------------------------------------------------------
for (dgp in dgps) {

  # Load and stack all jobs
  MC_all <- matrix(0, nrow = nrep * n_job, ncol = 60)
  for (job in 0:(n_job - 1)) {
    fmn <- paste0(address0, "/mc-dr.dgp-", dgp, ".n-", n, ".job-", job, ".RData")
    load(fmn)
    MC_all[(job * nrep + 1):((job + 1) * nrep), ] <- mc
  }

  #---------------------------------------------------------------------------
  # Summary statistics
  #---------------------------------------------------------------------------
  mean.mc        <- base::colMeans(MC_all, na.rm = TRUE)
  median.mc      <- base::apply(MC_all, 2, FUN = median, na.rm = TRUE)
  bias.mc        <- base::colMeans(MC_all - MC_all[, 1], na.rm = TRUE)
  median.bias.mc <- base::apply(MC_all - MC_all[, 1], 2, FUN = median, na.rm = TRUE)
  sd.mc          <- (base::colMeans(MC_all^2, na.rm = TRUE) -
                       base::colMeans(MC_all, na.rm = TRUE)^2)^0.5
  rmse.mc        <- base::colMeans((MC_all - MC_all[1])^2, na.rm = TRUE)^0.5
  mae.mc         <- base::colMeans(abs(MC_all - MC_all[1]), na.rm = TRUE)

  true_att      <- round(MC_all[1, 1], 3)
  sem_eff_bound <- round(mean.mc[35], 3)

  #---------------------------------------------------------------------------
  # Build panels
  #---------------------------------------------------------------------------
  panel_a <- build_est_rows(twfe_idx, mean.mc, median.bias.mc, rmse.mc)
  panel_b <- build_est_rows(dr_idx, mean.mc, median.bias.mc, rmse.mc)
  panel_c <- build_test_rows(test_idx, mean.mc)

  #---------------------------------------------------------------------------
  # Write CSV: header info + three panels separated by label rows
  #---------------------------------------------------------------------------
  out_file <- file.path(out_dir, paste0("Table ", dgp, ".csv"))

  # Determine column counts for alignment
  n_est_cols  <- ncol(panel_a)   # 7 columns
  n_test_cols <- ncol(panel_c)   # 5 columns
  max_cols    <- max(n_est_cols, n_test_cols)

  # Build header info rows
  header1 <- c(dgp_labels[as.character(dgp)],
               paste0("True ATT = ", true_att),
               paste0("Sem. Eff. Bound = ", sem_eff_bound),
               rep("", max_cols - 3))
  header2 <- c(paste0("n = ", n),
               paste0("MC reps = ", nrep * n_job),
               rep("", max_cols - 2))

  # Estimation column names
  est_header <- c("Method", "Avg. Bias", "Med. Bias", "RMSE", "Asy. Var.", "Cover.", "CIL")

  # Test column names (pad to max_cols)
  test_header <- c("Method", "Avg. Test Stats.", "Emp. Rej. Freq. (0.10)",
                    "Emp. Rej. Freq. (0.05)", "Emp. Rej. Freq. (0.01)",
                    rep("", max_cols - n_test_cols))

  # Panel label rows (pad to max_cols)
  panel_a_label <- c("Panel A: TWFE Estimators", rep("", max_cols - 1))
  panel_b_label <- c("Panel B: Nonparametric DR DiD Estimators", rep("", max_cols - 1))
  panel_c_label <- c("Panel C: Hausman-type Test", rep("", max_cols - 1))
  blank_row     <- rep("", max_cols)

  # Pad test rows to max_cols
  pad_test <- function(df) {
    m <- as.matrix(df)
    cbind(m, matrix("", nrow = nrow(m), ncol = max_cols - ncol(m)))
  }

  # Assemble full table as character matrix
  tbl <- rbind(
    header1,
    header2,
    blank_row,
    panel_a_label,
    est_header,
    as.matrix(panel_a),
    blank_row,
    panel_b_label,
    est_header,
    as.matrix(panel_b),
    blank_row,
    panel_c_label,
    test_header,
    pad_test(panel_c)
  )

  write.table(tbl, file = out_file, sep = ",", row.names = FALSE,
              col.names = FALSE, quote = TRUE)

  cat("DGP", dgp, "->", out_file, "\n")
}
