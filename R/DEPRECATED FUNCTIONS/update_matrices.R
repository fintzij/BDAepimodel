#' Update tpms, tpm products, emission matrices, and FB matrices before drawing
#' a new subject path.
#' 
#' @inheritParams draw_trajec
#'   
#' @return updated matrices within the epimodel environment.
#' 
#' @export

update_matrices <- function(epimodel, subject) {
          
          # update the tpms
          build_tpm_seqs(epimodel)

          # update the emission probability matrix
          build_emission_mat(epimodel)
          
          # construct the forward-backward matrices
          build_fb_mats(epimodel)
          
}