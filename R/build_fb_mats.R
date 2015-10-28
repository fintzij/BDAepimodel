#' Construct the forward-backward matrices in the forward pass of the FB
#' algorithm.
#' 
#' @param fb_mats array of forward-backward matrices to be updated
#' @param epimodel
#'   
#' @return list containing the FB matrices for the stochastic FB algorithm.
#' @export

build_fb_mats <- function(fb_mats, epimodel) {
          
          # calculate the emission probability-weighted initial distribution
          pi_0 <- normalize(epimodel$initdist * epimodel$emission_mat[,1, drop = FALSE])
          
          # compute the first fb matrix

          fb_mats[,,1] <- normalize((pi_0 %*% t(epimodel$emission_mat[, 2, drop = FALSE])) * epimodel$tpm_products[,,1])
          
          for(m in 2:length(epimodel$fb_mats)) {
                    
                    fb_mats[,,m] <- normalize((colSums(epimodel$fb_mats[,,m-1]) %*% t(epimodel$emission_mat[, m+1, drop = FALSE])) * epimodel$tpm_products[,,epimodel$.obs_time_inds[m]])
                    
          }
          
          return(fb_mats)
          
}