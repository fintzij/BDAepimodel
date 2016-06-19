#' Build a vector of probabilities according to which subjects should be
#' resampled.
#' 
#' @param epimodel
#'   
#' @return vector of probabilities according to which subject IDs are sampled
#' @export
subject_probabilities <- function(epimodel) {
          
          if(is.null(epimodel$sim_settings$preferential_sampling)) {
                    subj_probs <- normalize(rep(1, epimodel$popsize))
          } else {
                    pref_size   <- epimodel$sim_settings$preferential_sampling[[1]]
                    pref_weight <- epimodel$sim_settings$preferential_sampling[[2]]
                    
                    if(epimodel$sim_settings$preferential_sampling[[3]]) {
                              
                              # identify subjects with nonconstant paths
                              nonconstant_paths <- which(as.logical(match(seq_len(epimodel$popsize), epimodel$pop_mat[,"ID"], nomatch = FALSE)))
                              
                              # get subject IDs for the correct number of preferentially sampled subjects
                              if(length(nonconstant_paths) >= pref_size) {
                                        pref_subjects <- sort(sample(nonconstant_paths, pref_size))
                                        nonpref_subjects <- setdiff(seq_len(epimodel$popsize), pref_subjects)
                              } else {
                                        pref_subjects <- sort(c(nonconstant_paths, sample(setdiff(seq_len(epimodel$popsize), nonconstant_paths), pref_size - length(nonconstant_paths))))
                                        nonpref_subjects <- setdiff(seq_len(epimodel$popsize), pref_subjects)
                              }
                              
                              # assign probabilities
                              subj_probs <- rep(0, epimodel$popsize)
                              subj_probs[pref_subjects] <- normalize(rep(1, pref_size))* pref_weight
                              subj_probs[nonpref_subjects] <-  normalize(rep(1, epimodel$popsize - pref_size))* (1-pref_weight)
                              
                    } else {
                              subj_probs <- c(normalize(rep(1, pref_size))* pref_weight, normalize(rep(1, epimodel$popsize - pref_size))* (1-pref_weight))
                    }
          }
          
          return(subj_probs)
          
}