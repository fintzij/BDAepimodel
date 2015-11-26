#' Extracts a list of functions for computing flow rates.
#' 
#' @param states vector of compartment names
#' @param param_names vector of parameter names
#' @inheritParams simulate_epimodel
#'   
#' @return Returns a list of functions to compute rates. 
#' @export
#' 
#' @examples rates <- c("beta * I", "mu")
#' states <- c("S", "I", "R")
#' param_names <- c("beta", "mu")
#' extractRateFunctions(rates, states, param_names)      
#'        

extractRateFunctions <- function(rates, states, param_names) {
          
          for(t in 1:length(rates)){
                    
                    if(any(sapply(states, grepl, rates[t]))){
                              orig <- as.list(states[sapply(states, grepl, rates[t])])
                              repl <- as.list(paste("state(j,", which(sapply(states, grepl, rates[t])) + 2, ")", sep = ""))
                              for(s in 1:length(orig)){
                                        rates[t] <- gsub(orig[[s]], repl[[s]], rates[t])
                              }
                    }
                    
                    if(any(sapply(param_names, grepl, rates[t]))){
                              orig <- as.list(param_names[sapply(param_names, grepl, rates[t])])
                              repl <- as.list(paste("params[", which(sapply(param_names, grepl, rates[t])) - 1, "]", sep = ""))
                              
                              for(s in 1:length(orig)){
                                        rates[t] <- gsub(orig[[s]], repl[[s]], rates[t])
                              }
                    }
          }
          
          # initialize list of rate functions
          .rate_fcns <- vector("list", length = length(rates))
          
          # populate list with functions
          for(.s in seq_along(.rate_fcns)) {
                    
                    # capture text for the body of the function
                    .rate_fcns[[.s]] <- Rcpp::cppFunction(paste(paste0("Rcpp::NumericVector rate", .s),"(Rcpp::NumericMatrix& state, Rcpp::NumericVector& params) {
          Rcpp::IntegerVector stateDims = state.attr(\"dim\");
          Rcpp::NumericVector rate(stateDims[0]);
          for(int j=0; j < stateDims[0]; ++j) {
                    rate[j] = ",rates[.s],";
          }
          return rate;}"))
          } 
          
          return(.rate_fcns)
}