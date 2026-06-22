#-----------------------------------------------------------------------------
# This file gathers results from various job batches
# and creates the power curve plot (Figure 1).
#
# Note: dgp0 and delta0 are used as loop/parameter variables because the
# loaded RData files contain variables named `dgp` and `delta` that would
# overwrite them.
#-----------------------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)

#-----------------------------------------------------------------------------
# Parameters
#-----------------------------------------------------------------------------
n_job  <- 1
n      <- 1000
nrep   <- 500       # MC replications per job
dgp0   <- 3

delta_list <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
                0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)

# Paths
address0 <- "./simulation/Power Analysis (Figure 1)"  # repo-root-relative; run from the repository root
results_dir <- file.path(address0, "results")
funs_dir    <- file.path(address0, "funs")
fig_dir     <- file.path(address0, "figures")

##############################################################################
###############          LOAD PLOT DATA                           ############
##############################################################################

# Accumulators for test results
tests.ml.stb <- tests.ls.stb <- tests.ml.unstb <- tests.ls.unstb <- data.frame(
  delta = numeric(),
  avg_stat = numeric(),
  freq_10 = numeric(),
  freq_5 = numeric(),
  freq_1 = numeric(),
  stringsAsFactors = FALSE
)

for (delta0 in delta_list) {
  # Load job-level RData
  fmn <- paste0(results_dir, "/rlst_power.dgp-", dgp0,
                ".delta-", delta0, ".n-", n, ".job-0.RData")
  load(fmn)

  # Stack MC replications (single job)
  MC_all <- mc[1:nrep, ]

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

  #---------------------------------------------------------------------------
  # Summarize test results
  #---------------------------------------------------------------------------
  # Stabilized weights
  tests.ml.stb[nrow(tests.ml.stb) + 1, ] <- c(delta0,
                                                mean.mc[36],
                                                mean.mc[43],
                                                mean.mc[44],
                                                mean.mc[45])

  tests.ls.stb[nrow(tests.ls.stb) + 1, ] <- c(delta0,
                                                mean.mc[37],
                                                mean.mc[43],
                                                mean.mc[44],
                                                mean.mc[45])

  # Unstabilized weights
  tests.ml.unstb[nrow(tests.ml.unstb) + 1, ] <- c(delta0,
                                                    mean.mc[38],
                                                    mean.mc[46],
                                                    mean.mc[47],
                                                    mean.mc[48])

  tests.ls.unstb[nrow(tests.ls.unstb) + 1, ] <- c(delta0,
                                                    mean.mc[39],
                                                    mean.mc[49],
                                                    mean.mc[50],
                                                    mean.mc[51])
}


################################################################################
#############         MAKE POWER CURVE PLOTS                        ############
################################################################################
source(file.path(dirname(address0), "core", "dgps_drdid.R"))

kl_divg <- rep(0, length(delta_list))
i <- 0
for (delta in delta_list){
  i <- i + 1
  data  <- dgps_did(delta = delta, n = 1000000, xi_ps = 1)
  probs_p <- data$pi
  probs_q <- data$pp
  kl_divg[i] <- (mean(rowSums((log(probs_p) - log(probs_q)) * probs_p)))
}


df1 <- cbind(tests.ml.stb, kl_divg)
df2 <- cbind(tests.ls.stb, kl_divg)
df3 <- cbind(tests.ml.unstb, kl_divg)
df4 <- cbind(tests.ls.unstb, kl_divg)

# Add an identifier for each data frame
df1$Source <- "ML"
df2$Source <- "LS"


# Combine all data frames into a single data frame
combined_data <- rbind(
  melt(df1, id.vars = c("kl_divg", "Source"), measure.vars = c("freq_10", "freq_5", "freq_1"),
       variable.name = "Metric", value.name = "Value"),
  melt(df2, id.vars = c("kl_divg", "Source"), measure.vars = c("freq_10", "freq_5", "freq_1"),
       variable.name = "Metric", value.name = "Value")
)


# Rename metrics for clarity
combined_data$Metric <- recode(
  combined_data$Metric,
  freq_10 = "10%",
  freq_5 = "5%",
  freq_1 = "1%"
)
# Define a modern color palette for the metrics
custom_colors <- c("10%" = "#D81B60",
                   "5%" = "#1E88E5",
                   "1%" = "#FFC107") #004D40


# Custom line types for data sources
custom_linetypes <- c(
  "ML" = "solid",
  "LS" = "dotted"
)



p1 <- ggplot(combined_data, aes(x = kl_divg, y = Value, color = Metric, linetype = Source)) +
  geom_line(size = 1) +
  geom_point(shape = 15, size = 1.5) +
  labs(
    x = "Mean K-L Divergence",
    y = "Empirical Rejection Rate",
    color = "Nominal Test Size",
    linetype = "Data Source"
  ) +
  scale_color_manual(values = custom_colors) +
  scale_linetype_manual(values = custom_linetypes) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 15)
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 10)
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  guides(
    color = guide_legend(
      title = "Nominal Test Size",
      override.aes = list(linetype = "solid")
    ),
    linetype = guide_legend(
      title = "Cross-Validation Algorithm",
      override.aes = list(color = "black")
    )
  )

ggsave(file.path(fig_dir, "Figure 1.pdf"),
       plot = p1,
       width = 10, height = 5,
       bg = "transparent")
