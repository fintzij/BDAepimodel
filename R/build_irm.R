#' Construct an environment of rate matrices and another environment for the
#' corresponding eigen decompositions.
#'
#' Each rate matrix in the array has dimension \code{num_states x num_states}.
#'
#' @inheritParams simulate_epimodel
#'
#' @return environment with rate matrices
#' @export
#'
build_irm <- function(epimodel) {

        # compute the rates
        rates <- do.call(cbind, lapply(epimodel$rates, do.call, list(state = epimodel$irm_key_lookup, params = epimodel$params)))

        # update the rate matrices
        buildRateArray(irm_array = epimodel$irm, rates = rates, flow_inds = epimodel$flow_inds)

        # update the eigen decompositions
        if(is.null(epimodel$sim_settings$analytic_eigen)) {
                  buildEigenArray(real_eigenvals = epimodel$real_eigen_values,
                                  imag_eigenvals = epimodel$imag_eigen_values,
                                  eigenvecs      = epimodel$eigen_vectors, 
                                  inversevecs    = epimodel$inv_eigen_vectors, 
                                  irm_array      = epimodel$irm, 
                                  n_real_eigs    = epimodel$n_real_eigs)
        } else {
                  do.call(epimodel$sim_settings$analytic_eigen, 
                          args = list(real_eigenvals = epimodel$real_eigen_values,
                                      imag_eigenvals = epimodel$imag_eigen_values,
                                      eigenvecs      = epimodel$eigen_vectors, 
                                      inversevecs    = epimodel$inv_eigen_vectors, 
                                      irm_array      = epimodel$irm, 
                                      n_real_eigs    = epimodel$n_real_eigs))
        }
        
        return(epimodel)

}