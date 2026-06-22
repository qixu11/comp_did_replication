#-----------------------------------------------------------------------------
# This file gathers results from various job batches
# and creates a unified table using the merged results across all DGPs.
#-----------------------------------------------------------------------------
rm(list = ls())
n_job  <- 10
n      <- 1000
dgps   <- c(1, 3, 2)

address0 <- "./simulation/"
in_dir   <- paste0(address0, "Bandwidth Choice (Table 3)/results/") # EDIT to your local simulation path
out_dir  <- paste0(address0, "Bandwidth Choice (Table 3)/tables/")  # EDIT to your local simulation path

bw_methods   <- c("loocv",  "rcv", "plugin")
est_methods  <- c("nstnr_ml", "nstnr_ls", "stnr_ml", "stnr_ls")
test_methods <- c("stb_ml", "stb_ls")

# Accumulators for the unified tables
all_est  <- NULL
all_test <- NULL

#-----------------------------------------------------------------------------
# Loop over DGPs and bandwidth methods
#-----------------------------------------------------------------------------
for (dgp0 in dgps) {
  for (bw_method in bw_methods) {

    # Load and stack all jobs
    MC_all <- matrix(0, nrow = 100 * n_job, ncol = 60)
    for (job in 0:9) {
      fmn <- paste0(in_dir, bw_method, "/rlst_", bw_method,
                     ".dgp-", dgp0, ".n-1000.job-", job, ".RData")
      load(fmn)
      MC_all[(job * 100 + 1):((job + 1) * 100), ] <- mc[1:100, ]
    }

    #-------------------------------------------------------------------------
    # Summary statistics (identical to original logic)
    #-------------------------------------------------------------------------
    mean.mc       <- base::colMeans(MC_all, na.rm = TRUE)
    median.mc     <- base::apply(MC_all, 2, FUN = median, na.rm = TRUE)
    bias.mc       <- base::colMeans(MC_all - MC_all[, 1], na.rm = TRUE)
    median.bias.mc <- base::apply(MC_all - MC_all[, 1], 2, FUN = median, na.rm = TRUE)
    sd.mc         <- (base::colMeans(MC_all^2, na.rm = TRUE) -
                        base::colMeans(MC_all, na.rm = TRUE)^2)^0.5
    rmse.mc       <- base::colMeans((MC_all - MC_all[1])^2, na.rm = TRUE)^0.5
    mae.mc        <- base::colMeans(abs(MC_all - MC_all[1]), na.rm = TRUE)

    #-------------------------------------------------------------------------
    # Estimation results: 4 estimators per (dgp, bw_method)
    #-------------------------------------------------------------------------
    # Column indices into mc for each estimator:
    #   nstnr_ml: est=3, avar=4, cov=5, cil=6
    #   nstnr_ls: est=7, avar=8, cov=9, cil=10
    #   stnr_ml:  est=11, avar=12, cov=13, cil=14
    #   stnr_ls:  est=15, avar=16, cov=17, cil=18
    #   semipar eff bound: col 35
    est_idx <- list(
      nstnr_ml = c(3, 4, 5, 6),
      nstnr_ls = c(7, 8, 9, 10),
      stnr_ml  = c(11, 12, 13, 14),
      stnr_ls  = c(15, 16, 17, 18)
    )

    for (m in names(est_idx)) {
      idx <- est_idx[[m]]
      row <- data.frame(
        DGP            = dgp0,
        n              = n,
        Bandwidth      = bw_method,
        Method         = m,
        Av.Bias        = mean.mc[idx[1]] - mean.mc[1],
        Med.Bias       = median.bias.mc[idx[1]],
        RMSE           = rmse.mc[idx[1]],
        Asy.Var        = mean.mc[idx[2]],
        Coverage       = mean.mc[idx[3]],
        Length.CI      = mean.mc[idx[4]],
        Sem.Eff.Bound  = mean.mc[35],
        stringsAsFactors = FALSE
      )
      all_est <- rbind(all_est, row)
    }

    #-------------------------------------------------------------------------
    # Test results: 2 estimators per (dgp, bw_method)
    #-------------------------------------------------------------------------
    # stb_ml: stat=36, rej10=40, rej05=41, rej01=42
    # stb_ls: stat=37, rej10=43, rej05=44, rej01=45
    test_idx <- list(
      stb_ml = c(36, 40, 41, 42),
      stb_ls = c(37, 43, 44, 45)
    )

    for (m in names(test_idx)) {
      idx <- test_idx[[m]]
      row <- data.frame(
        DGP                  = dgp0,
        n                    = n,
        Bandwidth            = bw_method,
        Method               = m,
        Av.Stats             = mean.mc[idx[1]],
        Emp.Rej.Freq.0.10    = mean.mc[idx[2]],
        Emp.Rej.Freq.0.05    = mean.mc[idx[3]],
        Emp.Rej.Freq.0.01    = mean.mc[idx[4]],
        stringsAsFactors = FALSE
      )
      all_test <- rbind(all_test, row)
    }
  }
}

#-----------------------------------------------------------------------------
# Round and write unified tables
#-----------------------------------------------------------------------------
num_cols_est  <- c("Av.Bias", "Med.Bias", "RMSE", "Asy.Var",
                   "Coverage", "Length.CI", "Sem.Eff.Bound")
num_cols_test <- c("Av.Stats", "Emp.Rej.Freq.0.10",
                   "Emp.Rej.Freq.0.05", "Emp.Rej.Freq.0.01")

all_est[num_cols_est]   <- round(all_est[num_cols_est], 3)
all_test[num_cols_test] <- round(all_test[num_cols_test], 3)

write.csv(all_est,  file = file.path(out_dir, paste0("bw.estimator.summary.n-", n, ".csv")),
          row.names = FALSE)
write.csv(all_test, file = file.path(out_dir, paste0("bw.test.summary.n-", n, ".csv")),
          row.names = FALSE)

cat("Estimation table:", nrow(all_est), "rows written.\n")
cat("Test table:      ", nrow(all_test), "rows written.\n")


#-----------------------------------------------------------------------------
# Build a combined table: estimation + test columns in one table
#-----------------------------------------------------------------------------
rows_est <- c(1, 5, 9, 13, 17, 21, 25, 29, 33) 
rows_test <- c(1, 3, 5, 7, 9, 11, 13, 15, 17) 
# 1) Keep only the columns you need from each table
est_keep <- all_est[rows_est, c("DGP", "n", "Bandwidth", "Method",
                        "Av.Bias", "Med.Bias", "RMSE", "Coverage", "Length.CI")]

test_keep <- all_test[rows_test, c("DGP", "n", "Bandwidth", "Method",
                          "Av.Stats",
                          "Emp.Rej.Freq.0.10", "Emp.Rej.Freq.0.05", "Emp.Rej.Freq.0.01")]


test_keep2 <- test_keep[, c("DGP","n","Bandwidth","Av.Stats",
                            "Emp.Rej.Freq.0.10","Emp.Rej.Freq.0.05","Emp.Rej.Freq.0.01")]

combined <- merge(est_keep, test_keep2, by = c("DGP","n","Bandwidth"), all.x = TRUE)

# 3) Rename columns to your preferred labels
names(combined)[names(combined) == "Av.Bias"]             <- "Avg. Bias"
names(combined)[names(combined) == "Med.Bias"]            <- "Med. Bias"
names(combined)[names(combined) == "Coverage"]            <- "Cover."
names(combined)[names(combined) == "Length.CI"]           <- "CIL"
names(combined)[names(combined) == "Av.Stats"]            <- "Avg. Test Stats."
names(combined)[names(combined) == "Emp.Rej.Freq.0.10"]    <- "Emp. Rej. Freq. (0.10)"
names(combined)[names(combined) == "Emp.Rej.Freq.0.05"]    <- "Emp. Rej. Freq. (0.05)"
names(combined)[names(combined) == "Emp.Rej.Freq.0.01"]    <- "Emp. Rej. Freq. (0.01)"

# 4) Keep only the requested columns (plus identifiers so it’s usable)
combined_out <- combined[, c("DGP","n","Bandwidth","Method",
                             "Avg. Bias","Med. Bias","RMSE","Cover.","CIL",
                             "Avg. Test Stats.",
                             "Emp. Rej. Freq. (0.10)",
                             "Emp. Rej. Freq. (0.05)",
                             "Emp. Rej. Freq. (0.01)")]

dgp_labels <- c(
  "Non-Stationary Covariate Distribution (delta^alt = 0)",
  "Non-Stationary Covariate Distribution (delta^alt = 0.5)",
  "Stationary Covariate Distribution (delta^alt = 1)"
)

# Map numeric DGP to labels
combined_out$DGP_Label <- dgp_labels[match(combined_out$DGP, c(1,3,2))]
combined_out$DGP_Label <- factor(
  combined_out$DGP_Label,
  levels = dgp_labels
)
combined_out$Bandwidth <- factor(combined_out$Bandwidth, levels = bw_methods)


combined_out <- combined_out[order(combined_out$DGP_Label, combined_out$Bandwidth), ]

# 6) Write it
write.csv(combined_out,
          file = file.path(out_dir, paste0("Table 3.csv")),
          row.names = FALSE)
