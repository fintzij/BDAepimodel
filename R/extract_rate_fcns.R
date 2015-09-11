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
#' extract_rate_fcns(rates, states, param_names)      
#'        

extract_rate_fcns <- function(rates, states, param_names) {
          
          require(pryr, quietly = TRUE)
          
          for(t in 1:length(rates)){
                    
                    rates[t] <- paste("unname(", rates[t], ")", sep = "")
                    
                    if(any(sapply(states, grepl, rates[t]))){
                              orig <- as.list(states[sapply(states, grepl, rates[t])])
                              repl <- as.list(paste("state[,\"", states[sapply(states, grepl, rates[t])], "\"]", sep = ""))
                              for(s in 1:length(orig)){
                                        rates[t] <- gsub(orig[[s]], repl[[s]], rates[t])
                              }
                    }
                    
                    if(any(sapply(param_names, grepl, rates[t]))){
                              orig <- as.list(param_names[sapply(param_names, grepl, rates[t])])
                              repl <- as.list(paste("params[\"", param_names[sapply(param_names, grepl, rates[t])], "\"]", sep = ""))
                              
                              for(s in 1:length(orig)){
                                        rates[t] <- gsub(orig[[s]], repl[[s]], rates[t])
                              }
                    }
          }

          # convert rates to list of unevaluated expressions
          .rates <- eval(parse(text = paste("alist(",paste(c(rates), collapse = ", "), ")", sep = "")))
          
          # initialize list of rate functions
          .rate_fcns <- vector("list", length = length(rates))
          
          # populate list with functions
          for(.s in seq_along(.rate_fcns)) {
                    
                    # generate a list of arguments
                    .arg_list <- eval(parse(text = paste("alist(",paste(c("state", "params","...="), collapse = "=,"),")",sep="")))
                    
                    # capture text for the body of the function
                    .rate_fcns[[.s]] <- make_function(.arg_list, .rates[[.s]])
          }
          
          return(.rate_fcns)
}