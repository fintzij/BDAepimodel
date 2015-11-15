#' Template for function for transformint to the estimation scale
#' 
#' Template for writing a MH transition kernel, example based on SIR model, 
#' where proposals for the infectivity and recovery parameters are made in terms
#' of log(R_0), and log(mu). See package vignette for additional details. The
#' function for transforming from the estimation scale back to the model scale
#' should be similarly specified.
#' 
#' @param epimodel
#'   
#' @return template for transformation functions.
#' @export
#' 
to_est_scale_template <- function(params, epimodel) {
          
          ############################################################################
          ### COPY AND PASTE THIS TEMPLATE INTO A NEW SCRIPT AND EDIT ACCORDINGLY. ###
          ### DO NOT EDIT IN PLACE!!!!                                             ###
          ############################################################################
          
          # get parameters
          params_est          <- params
          
          # transform the parameters by indexing into the appropriate elements
          # by name. Any object within an epimodel list may be accessed.
          params_est["beta"]  <- log(params["beta"] * epimodel$popsize/params["mu"])
          
          params_est["mu"]    <- log(params["mu"])
          
          # return a vector of parameters on the estimation scale. 
          return(params_est)
}
