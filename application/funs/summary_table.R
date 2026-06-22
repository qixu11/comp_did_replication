#-----------------------------------------------------------------------------------------------
# Create Table 4 from saved estimation results.
#
# Reads:  results/result_sequeira_drdid.RData (DR DiD estimates)
#         results/result_sequeira_fe.RData (TWFE estimates)
# Writes: tables/Table 4.xlsx    (two sheets: "est" and "test")
#
# Sheet "est": estimation results for TWFE_1, TWFE_2, DR_stnr, DR_nstnr
#   Each block shows: ATT, asymptotic SE in (), cluster-robust SE in [].
#
# Sheet "test": Hausman test for covariate stationarity
#   Asymptotic and cluster-bootstrap p-values and rejection decisions.
#-----------------------------------------------------------------------------------------------

# Set the base directory to the application root.
# Users should modify this path to point to their local copy.
base_dir <- "~/Desktop/DID CC Claude/application"

source(file.path(base_dir, "funs/brackets.R"))

load(file.path(base_dir, "results/result_sequeira_drdid.RData"))
load(file.path(base_dir, "results/result_sequeira_twfe.RData"))

twfe.1.summary <- twfe.2.summary <- dml.summary <- dr.nstnr.summary <- dr.stnr.summary <- stnr.test.summary  <- out_table<-list()

for (i in 1:4){
  dr.nstnr.summary[[i]] <- rbind(round(sequeira_res[[i]]$att_nstnr,3),
                                 parens(round(sequeira_res[[i]]$se_nstnr,3)),
                                 bracket(round(sequeira_res[[i]]$sse_nstnr,3))

  ) 
  
  dr.stnr.summary[[i]] <- rbind(round(sequeira_res[[i]]$att_stnr,3),
                                parens(round(sequeira_res[[i]]$se_stnr,3)),
                                bracket(round(sequeira_res[[i]]$sse_stnr,3))

  ) 
  
  twfe.1.summary[[i]] <- rbind(round(sequeira_fe[[i]]$att_twfe1,3),
                               parens(round(sequeira_fe[[i]]$se_twfe1,3)),
                               bracket(round(sequeira_fe[[i]]$sse_twfe1,3))

  ) 
  
  twfe.2.summary[[i]] <- rbind(round(sequeira_fe[[i]]$att_twfe2,3),
                               parens(round(sequeira_fe[[i]]$se_twfe2,3)),
                               bracket(round(sequeira_fe[[i]]$sse_twfe2,3))

  ) 
  

  
}




header <-  c("Estimator", " ", "Prob(bribe)",
  "log(bribe)",
  "log(bribe/shpt.val.)",
  "log(bribe/shpt.tonn.)"
  )

coln_label <- cbind(c("TWFE_1", " ",  " ", 
                      "TWFE_2", " ",  " ",   
                      "DR_stnr", " ", " ", 
                      "DR_nstnr", " ", " "),  
                    c("att", "se.", "cls. se.",
                      "att", "se.", "cls. se.",
                      "att", "se.", "cls. se.",
                      "att", "se.", "cls. se."))





  summary_est_table <- rbind(cbind( twfe.1.summary[[1]], twfe.1.summary[[2]], twfe.1.summary[[3]], twfe.1.summary[[4]]),
                       cbind(twfe.2.summary[[1]], twfe.2.summary[[2]], twfe.2.summary[[3]], twfe.2.summary[[4]]),
                      cbind( dr.stnr.summary[[1]], dr.stnr.summary[[2]], dr.stnr.summary[[3]], dr.stnr.summary[[4]]),
                      cbind(dr.nstnr.summary[[1]], dr.nstnr.summary[[2]], dr.nstnr.summary[[3]], dr.nstnr.summary[[4]])
                          )
  out_table[[1]] <- rbind( header, cbind(coln_label, summary_est_table))
  
  
  
  
# summarize test results
  for (i in 1:4){
    stnr.test.summary[[i]] <- rbind(round(1-pchisq(sequeira_res[[i]]$wald_stat, df= 1),3),
                                    ifelse(sequeira_res[[i]]$stnr_test, "Yes","No"),
                                    round(sequeira_res[[i]]$p_value,3),
                                    ifelse(sequeira_res[[i]]$stnr_btest, "Yes","No")
    )
  }
  
  
  header_test <-  c("Estimator", 
                     "Prob(bribe)",
                    "log(bribe)",
                    "log(bribe/shpt.val.)",
                    "log(bribe/shpt.tonn.)"
  )                 
  coln_label_test <-  c("p_value", "H0_rejected?",
                        "cls_p_value", "cls_H0_rejected?")
  
  
  summary_test_table  <- cbind(stnr.test.summary[[1]],stnr.test.summary[[2]],
                               stnr.test.summary[[3]],stnr.test.summary[[4]]
                         )
    
    
  
  out_table[[2]] <- rbind( header_test, cbind(coln_label_test, summary_test_table))
  
  
# Save the output table  
library(openxlsx)
out_file <- file.path(base_dir, "tables/Table 4.xlsx")
out <- createWorkbook()
sname <- c("est", "test")

for (i in 1:2) {
  addWorksheet(out, sname[i])
  writeData(out, sheet = sname[i], colNames = FALSE, x = out_table[[i]])
}
saveWorkbook(out, out_file, overwrite = TRUE)

