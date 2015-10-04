#' Sample a subject-level trajectory.
#' 
#' Proceeds in four steps: 1) Sample disease status at observation times using
#' the stochastic forward-backward algorithm. 2) Sample the infection status at 
#' transition times of other subjects. 3) Sample exact times of transition 
#' within inter-event intervals where needed. 4) Accept/Reject the proposed 
#' trajectory via Metropolis-Hastings.
#' 
#' @inheritParams get_tpms_to_build
#'   
#' @return updated configuration matrix and observation matrix objects in the
#'   epimodel environment.
#' @export

draw_trajec <- function(epimodel, subject) {
          
          # generate subject ID
          .subj_ID <- paste0(".X", subject)
          
          # set variable for where to start storing additional state changes
          epimodel$.subj_row_ind <- epimodel$.ind_final_config + 1
          
          # sample status at observation times
          sample_at_obs_times(epimodel, subject = subject, subj_ID = .subj_ID)
          
          # sample status at event times
          sample_at_event_times(epimodel, subject = subject, subj_ID = .subj_ID)
          
          # sample paths in inter-event intervals
          sample_path(epimodel, subject = subject, subj_ID = .subj_ID)

          # insert trajectory back into the configuration matrix 
          insert_trajectory(epimodel, subject = subject, subj_ID = .subj_ID, reinsertion = FALSE)
          
          # accept or reject the proposed path via metropolis-hastings
          MH_accept_reject(epimodel, subject = subject, subj_ID = .subj_ID)
}