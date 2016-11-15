#' Append missing objects to the irm arrays and eigen decomposition objects.
#' Update irm_key_lookup matrix.
#'
#' @inheritParams simulate_epimodel
#'
#' @return epimodel with complete irm and eigen decomposition arrays.
#' @export

append_missing <- function(epimodel) {

        # identify missing keys
        missing_keys <- which(epimodel$keys == 0)

        # get the required configurations from the configuration matrix
        missing_configs <- unique(epimodel$pop_mat[missing_keys, epimodel$index_states, drop = FALSE])

        # compute the rates for the missing configurations
        rates <- do.call(cbind, lapply(epimodel$rates, do.call, list(state = missing_configs, params = epimodel$params)))

        # new keys
        new_keys <- nrow(epimodel$irm_key_lookup) + seq_len(nrow(missing_configs))

        # update the lookup matrix
        epimodel$irm_key_lookup <- rbind(epimodel$irm_key_lookup, cbind(new_keys, missing_configs))

        # insert the missing keys
        for(r in 1:nrow(missing_configs)) {

              # add the key to the vector of keys
                  key_inds <- apply(epimodel$pop_mat[missing_keys, epimodel$index_states, drop = FALSE], 1, 
                                    function(x) all(x == missing_configs[r, , drop = FALSE]))
                  epimodel$keys[missing_keys[key_inds]] <- new_keys[r]
        }

        # create empty arrays for irm, eigen vectors and inverse vectors, and matrix for eigen values
        new_rate_array  <- array(0.0, dim = c(length(epimodel$states), length(epimodel$states), length(new_keys)))
        new_eigen_array <- array(0.0, dim = c(length(epimodel$states), length(epimodel$states), length(new_keys)))
        new_inv_array   <- array(0.0, dim = c(length(epimodel$states), length(epimodel$states), length(new_keys)))
        new_real_values <- matrix(0.0, nrow = length(epimodel$states), ncol = length(new_keys))
        new_imag_values <- matrix(0.0, nrow = length(epimodel$states), ncol = length(new_keys))
        n_real_eigs     <- vector(mode = "integer", length = length(new_keys))

        # fill arrays and matrix
        buildRateArray(new_rate_array, rates, epimodel$flow_inds)
        
        if(is.null(epimodel$sim_settings$analytic_eigen)) {
                  
                  buildEigenArray(real_eigenvals = new_real_values,
                                  imag_eigenvals = new_imag_values,
                                  eigenvecs    = new_eigen_array,
                                  inversevecs  = new_inv_array,
                                  irm_array    = new_rate_array,
                                  n_real_eigs  = n_real_eigs
                  )
                  
        } else if(epimodel$sim_settings$analytic_eigen == "SIR") {
                  
                  buildEigenArray_SIR(real_eigenvals = new_real_values,
                                      imag_eigenvals = new_imag_values,
                                      eigenvecs      = new_eigen_array,
                                      inversevecs    = new_inv_array,
                                      irm_array      = new_rate_array,
                                      n_real_eigs    = n_real_eigs, 
                                      initial_calc   = TRUE
                  )
                  
        } else if(epimodel$sim_settings$analytic_eigen == "SEIR") {
                  
                  buildEigenArray_SEIR(real_eigenvals = new_real_values,
                                       imag_eigenvals = new_imag_values,
                                       eigenvecs      = new_eigen_array,
                                       inversevecs    = new_inv_array,
                                       irm_array      = new_rate_array,
                                       n_real_eigs    = n_real_eigs, 
                                       initial_calc   = TRUE
                  )
                  
        } 
        
        # insert new objects into existing arrays and matrix
        epimodel$irm               <- joinCubes(epimodel$irm, new_rate_array)
        epimodel$eigen_vectors     <- joinCubes(epimodel$eigen_vectors, new_eigen_array)
        epimodel$inv_eigen_vectors <- joinCubes(epimodel$inv_eigen_vectors, new_inv_array)
        epimodel$real_eigen_values <- cbind(epimodel$real_eigen_values, new_real_values)
        epimodel$imag_eigen_values <- cbind(epimodel$imag_eigen_values, new_imag_values)
        epimodel$n_real_eigs       <- c(epimodel$n_real_eigs, n_real_eigs)

        return(epimodel)
}