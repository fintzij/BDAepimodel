#' Instatiate missing IRMs, and update TPMS, TPM products, and FB matrices
#' before drawing a new subject path.
#' 
#' @inheritParams get_tpms_to_build
#' @param emission TRUE/FALSE for whether to construct the emission matrices
#' @param fb TRUE/FALSE for whether to construct the FB matrices
#'   
#' @return updated matrices within the epimodel environment.

update_matrices <- function(epimodel, subject, irm = TRUE, emission = TRUE, fb = TRUE) {
          
          # check to see if any additional irms are needed.
          # if so, check_irm will instatiate the required
          # matrices and their eigen decompositions
          if(irm) check_irm(epimodel)
          
          # get indices of tpms that need to be rebuilt
          epimodel$.tpms_to_build <- get_tpms_to_build(epimodel, subject = .subjects[j])
          
          # update the transition probability matrix sequences
          build_tpm_seqs(epimodel)
          
          # update the emission probability matrix
          if(emission) build_emission_mat(epimodel)
          
          # construct the forward-backward matrices
          if(fb) build_fb_mats(epimodel)
          
}