#' Initialize an arrays for storing matrices of rates, eigenvectors, and inverse matrices for eigenvectors, along with a matrix of eigenvalues.
#'
#' @param epimodel
#'
#' @return epimodel object with initialized objects
#' @export
#'
init_irm <- function(epimodel) {

        # initialize lookup matrix
        epimodel$irm_key_lookup <- as.matrix(unique(epimodel$pop_mat[1:epimodel$ind_final_config, epimodel$index_states, drop = FALSE]))
        epimodel$irm_key_lookup <- cbind(1:nrow(epimodel$irm_key_lookup), epimodel$irm_key_lookup)

        colnames(epimodel$irm_key_lookup)[1] <- "key"

        # instatiate the irm and eigen arrays
        n_states <- epimodel$num_states
        n_irms   <- nrow(epimodel$irm_key_lookup)
        
        epimodel$irm               <- array(0.0, dim = c(n_states, n_states, n_irms))
        epimodel$eigen_vectors     <- array(0.0, dim = c(n_states, n_states, n_irms))
        epimodel$inv_eigen_vectors <- array(0.0, dim = c(n_states, n_states, n_irms))
        epimodel$real_eigen_values <- matrix(0.0, nrow = n_states, ncol = n_irms)
        epimodel$imag_eigen_values <- matrix(0.0, nrow = n_states, ncol = n_irms)
        epimodel$n_real_eigs       <- vector(mode = "integer", length = n_irms)
        
        return(epimodel)
}