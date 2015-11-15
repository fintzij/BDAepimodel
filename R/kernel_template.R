#' Generates a template script for specifying a transition kernel or functions 
#' to and from the estimation scale.
#' 
#' @param type either "mh" for a Metropolis-Hastings kernel, "gibbs" for a Gibbs
#'   kernel, or "est_scale" for a transformation function from the model scale
#'   to the parameter estimation scale.
#'   
#' @return template
#' @export
#' 
templates <- function(type) {
          
          edit(eval(parse(text = paste0(type, "_template"))))
}