#' logit function
#' 
#' logit(x) = log(x/(1-x))
#'
#' @param x 
#'
#' @return logit of x
#'
logit <- function(x) {
          log(x/(1-x))
}

#' expit function
#' 
#' expit(x) = 1/(1+exp(x))
#'
#' @param x 
#'
#' @return expit of x
#' 
expit <- function(x) {
          1/(1 + exp(x))
}