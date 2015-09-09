#' Update the observation bookkeeping matrix with a new trajectory.
#' 
#' @inheritParams simulate_epimodel
#' @param direction update observation matrix to reflect the addition or removal
#'   of a subject-level trajectory. Valid arguments are \code{add} or 
#'   \code{remove}.
#' @param trajec matrix of configurations for a single subject, where each row
#'   is a row from \code{epimodel$config_mat} for that subject.
#'   
#' @return updated \code{obs_mat} matrix
#' 

update_obs_mat <- function(epimodel, direction, trajec) {
          
}