#' Updates the matrix of emission probabilities
#'
#' @param epimodel 
#'
#' @return updated matrix of emission probabilities
#' @export

build_emission_mat <- function(epimodel) {
          
          # calculate the emission probabilities for the unmeasured compartments
          epimodel$.emission_mat[epimodel$unmeasured_vars,]  <- matrix(epimodel$d_meas_process(state = epimodel$obs_mat, meas_vars = epimodel$meas_vars, params = epimodel$params, log = FALSE), nrow = epimodel$num_unmeasured, ncol = epimodel$nobs, byrow = TRUE)
          
          # for each of the measured compartments, calculate the emission probabilities
          for(m in 1:epimodel$num_measured) {
                    
                    .add_mat <- matrix(0, nrow = epimodel$nobs, ncol = 1 + 2*epimodel$num_measured)
                    .add_mat[,1 + 2*m] <- 1
                    
                    epimodel$.emission_mat[epimodel$meas_vars[m],] <- epimodel$d_meas_process(state = epimodel$obs_mat + .add_mat, meas_vars = epimodel$meas_vars, params = epimodel$params, log = FALSE)
                    
          }
          
}