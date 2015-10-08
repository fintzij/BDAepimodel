#' Construct a transition probability matrix from
#' 
#' If (t1 - t0) is small, there could be numerical underflow leading to negative
#' entries in the tpm. If so, the negative values are set to 0.
#' 
#' @param values vector of eigenvalues
#' @param vectors matrix of eigenvectors
#' @param inv_vectors inverse of matrix of eigenvectors
#' @param t0 left endpoint of the interval
#' @param t1 right endpoint of the interval
#'   
#' @return transition probability matrix
#' @export

build_tpm <- function(values, vectors, inv_vectors, t0, t1) {
          
          mat <- vectors%*%(.Internal(diag(exp(values*(t1-t0)), length(values), length(values)))%*%inv_vectors)
          
          # correct floating point errors
          if(any(mat < 0)){
                    mat[mat < 0] <- 0
          }
          
          return(mat)
}