#' Instatiate missing IRMs, and update TPMS, TPM products, and FB matrices
#' before drawing a new subject path.
#' 
#' @inheritParams get_tpms_to_build
#'   
#' @return updated matrices within the epimodel environment.

update_matrices <- function(epimodel, subject) {
          
          # check to see if any additional irms are needed.
          # if so, check_irm will instatiate the required
          # matrices and their eigen decompositions
          check_irm(epimodel)
          
          # get indices of tpms that need to be rebuilt
          epimodel$.tpms_to_build <- get_tpms_to_build(epimodel, subject = .subjects[j])
          
          # update the transition probability matrix sequences
          build_tpm_seqs(epimodel)
          
          # update the emission probability matrix
          build_emission_mat(epimodel)
          
          # construct the forward-backward matrices
          build_fb_mats(epimodel)
          
}