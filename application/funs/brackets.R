#' Formatting helpers for table output.
#'
#' bracket_ci: format confidence interval as "[lci, uci]"
#' bracket:    format scalar as "[x]" (used for cluster-robust SE)
#' parens:     format scalar as "(x)" (used for asymptotic SE)

bracket_ci <- function(lci, uci, digit = 3){
  lci <- round(lci, digit)
  uci <- round(uci, digit)
  
  ci_bracket = rep(NA, length(lci))
  for (i in 1:length(lci)){
    ci_bracket[i] = paste(c("[", lci[i], ",", uci[i], "]"), collapse="")
  }
  return(ci_bracket)
}


bracket <-function(x, digit = 3){
  x <- round(x, digit)
  bracket_x = paste(c("[", x,  "]"), collapse="")
  return(bracket_x)
}


parens <-function(x, digit = 3){
  x <- round(x, digit)
  parens_x = paste(c("(", x,  ")"), collapse="")
  return(parens_x)
}