#' Compute the eigen decomposition for a rate matrix, and either instatiate it 
#' or update an existing decomposition.
#' 
#' @param irm rate matrix
#' @inheritParams build_rate_mat
#'   
#' @return list containing the unordered eigenvalues, the matrix of
#'   eigenvectors, and the inverse of the eigenvector matrix
#' @export

irm_decomp <- function(irm, env = NULL, key = NULL, update = TRUE) {
          
          if(is.null(key) || update == FALSE) {
                    
                    # compute the eigen decomposition
                    .irm_decomp                   <- .Internal(La_rg(irm, FALSE))
                    .irm_decomp$inv_vectors       <- solve(.irm_decomp$vectors, diag(nrow(irm)))
                    
                    return(.irm_decomp)
                    
          } else {
                    env[[key]][1:2]               <- .Internal(La_rg(irm, FALSE))
                    env[[key]][[3]]               <- solve(env[[key]]$vectors, diag(nrow(irm)))
                    
          }
}