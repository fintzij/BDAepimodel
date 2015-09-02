#' Computes the lumped rates for simulation purposes.
#' 
#' @inheritParams init_epimodel
#' @param status_vec vector of subject states at the current time.
#' @param time current time used for indexing into the pop_mat matrix
#'   
#' @return vector of lumped rates with number of elements equal to the number of
#'   possible transitions.
#'   
# lump_rates <- function(epimodel, status_vec, time) {
#           
#           lumped_rates <- rep(0, ncol(epimodel$flow))
#           
#           for(r in 1:length(lumped_rates)){
#                     lumped_rates[r] <- 1
#           }
# }