#' Build a vector of probabilities according to which subjects should be
#' resampled.
#' 
#' @param epimodel
#'   
#' @return vector of probabilities according to which subject IDs are sampled
#' @export
subject_probabilities <- function(epimodel) {
          
          preferential_sampling <- !is.null(epimodel$sim_settings$preferential_sampling)
          
          if(!preferential_sampling) {
                    subj_probs <- normalize(rep(1, epimodel$popsize))
                    pref_size  <- NULL
                    pref_prob  <- NULL
                    pref_IDs   <- NULL
                    nonpref_IDs <- NULL
          } else {
                    if(epimodel$sim_settings$preferential_sampling[1] > epimodel$popsize) {
                              stop("The size of the preferentially sampled group must not be larger than the population size.")
                    }
                    
                    subj_probs  <- NULL
                    pref_size   <- epimodel$sim_settings$preferential_sampling[1] # size of the preferentially sampled group
                    pref_prob   <- epimodel$sim_settings$preferential_sampling[2] # preferential sampling probability
                    pref_IDs    <- setdiff(unique(epimodel$pop_mat[,"ID"]), c(0, NA)) # preferential sampling IDs
                    nonpref_IDs <- setdiff(seq_len(epimodel$popsize), pref_IDs)       # non-preferential sampling IDs
                    
                    # ensure the preferential group is the right size
                    if(length(pref_IDs) < pref_size) {
                              # sample IDs to add to the pref set
                              if(length(nonpref_IDs) != 1) {
                                        add_IDs = sample(nonpref_IDs, pref_size - length(pref_IDs)) 
                              } else {
                                        add_IDs = nonpref_IDs
                              }
                              pref_IDs    = c(pref_IDs, add_IDs)          # add IDs to the preferential set
                              nonpref_IDs = setdiff(nonpref_IDs, add_IDs) # remove IDs from the non-preferential set
                              
                    } else if(length(pref_IDs) > pref_size) {
                              # sample IDs to remove from the pref set
                              if(length(pref_IDs) != 1) {
                                        rem_IDs = sample(pref_IDs, length(pref_IDs) - pref_size) 
                              } else {
                                        rem_IDs = pref_IDs
                              }
                              pref_IDs    = setdiff(pref_IDs, rem_IDs) # remove IDs from the preferential set
                              nonpref_IDs = c(nonpref_IDs, rem_IDs)    # add IDs to the non-preferential set
                    }
          }
          
          samp_probs <- list(preferential_sampling = preferential_sampling,
                             pref_size             = pref_size,
                             pref_prob             = pref_prob,
                             subj_probs            = subj_probs,
                             pref_IDs              = pref_IDs,
                             nonpref_IDs           = nonpref_IDs)
          
          return(samp_probs)
}