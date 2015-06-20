#' Extracts a list of functions for computing flow rates.
#' 
#' @param compartments vector of compartment names
#' @param params vector of parameter names
#' @inheritParams simulate_epimodel
#'   
#' @return Exports list of rate functions to the global environment. 
#' @export
#' 
extract_rate_fcns <- function(rates, compartments, params) {
          
          # initialize list of rates
          .rate_fcns <- vector("list", length = length(rates))
          names(.rate_fcns) <- paste("rate",1:length(rates), sep = "")
          
          # populate list with functions
          for(.s in seq_along(.rate_fcns)) {
                    
                    # identify the compartments and parameters for the rate
                    .args_s <- c(compartments[sapply(compartments, grepl, rates[.s])], 
                                 params[sapply(params, grepl, rates[.s])])
                    
                    # generate a list of arguments
                    .arg_list <- vector("list", length = length(.args_s))
                    names(.arg_list) <- .args_s
                    
                    .rate_fcns[[.s]] <- make_function(args = as.pairlist(.arg_list), body = quote(eval(parse(text = rates[.s]))))
          }
          
          return(.rate_fcns)
}


