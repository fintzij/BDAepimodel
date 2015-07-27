#' Extracts a list of functions for computing flow rates.
#' 
#' @param compartments vector of compartment names
#' @param param_names vector of parameter names
#' @inheritParams simulate_epimodel
#'   
#' @return Returns a list of functions to compute rates. 
#' @export
#' 
#' @examples rates <- c("beta * S * I", "mu * I")
#' compartments <- c("S", "I", "R")
#' param_names <- c("beta", "mu")
#' extract_rate_fcns(rates, compartments, param_names)      
#'        

extract_rate_fcns <- function(rates, compartments, param_names) {
          
          # convert rates to list of unevaluated expressions
          .rates <- eval(parse(text = paste("alist(",paste(c(rates), collapse = ", "), ")", sep = "")))
          
          # initialize list of rate functions
          .rate_fcns <- vector("list", length = length(rates))
          
          # populate list with functions
          for(.s in seq_along(.rate_fcns)) {
                    
                    # identify the compartments and parameters for the rate
                    .args_s <- c(compartments[sapply(compartments, grepl, rates[.s])], 
                                 param_names[sapply(param_names, grepl, rates[.s])])
                    
                    # generate a list of arguments
                    .arg_list <- eval(parse(text = paste("alist(",paste(c(.args_s,"...="), collapse = "=,"),")",sep="")))
                    
                    # capture text for the body of the function
                    .rate_fcns[[.s]] <- make_function(.arg_list, .rates[[.s]])
          }

          
          return(.rate_fcns)
}
