#' Initializes a list of bookkeeping matrices for the observations, population 
#' level trajectory, and subject level trajectory.
#' 
#' @param dat logical indicator for whether to include the data matrix in the 
#'   bookkeeping list.
#' @param pop logical indicator for whether to include the population level 
#'   trajectory matrix in the bookkeeping list.
#' @param subj logical indicator for whether to include the subject level
#'   trajectory matrix in the bookkeeping list.
#' @inheritParams simulate_epimodel
#'   
#' @return list of matrices for data, population trajectory, and subject level trajectory
#' @export
#' 
init_bookkeeping <- function(dat = TRUE, pop = TRUE, subj = TRUE, obstimes, meas_vars, init_state, tmax) {
          
          # initialize list object
          epimod_mats <- vector("list", length = sum(c(dat, pop, subj)))
          names(epimod_mats) <- c("data", "population", "subjects")[c(dat, pop, subj)]
          
          # if data is TRUE, set up data matrix
          if(dat) epimod_mats$data <- init_data_mat(obstimes = obstimes, meas_vars = meas_vars)
          
          # if population is TRUE, set up population level trajectory matrix
          if(pop) epimod_mats$population <- init_pop_mat(init_state = init_state, tmax = tmax)
          
          # if subj_mat is TRUE (default), set up individual level trajectory object
          if(subj) epimod_mats$subjects <- init_subj_mat(init_state = init_state, tmax = tmax)
          
          # return list of matrices
          return(epimod_mats)
}