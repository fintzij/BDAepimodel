#' Calculate the likelihood for a subject-level trajectory based on a
#' time-inhomogeneous CTMC.
#'
#' @inheritParams calc_pop_likelihood
#' @param subj_ID string for subject ID
#'
#' @return likelihood or log-likelihood
#' @export

calc_subj_likelihood <- function(epimodel, subject, log = TRUE) {

        # get the appropriate indices to exclude the intermediate observation time rows
        .left_endpoints     <- c(1, c(1:(epimodel$ind_final_config))[-epimodel$obs_time_inds])
        .right_endpoints    <- c(.left_endpoints[-1], epimodel$ind_final_config)

          # get the subject_path
          .subj_path          <- epimodel$config_mat[c(1, .right_endpoints), subject]
          .path_length        <- length(.subj_path)

          # compute the time diffs
          .time_diffs <- epimodel$pop_mat[.right_endpoints, "time"] - epimodel$pop_mat[.left_endpoints, "time"]

          # get irm keys
          .irm_keys <- retrieveKeys(.left_endpoints, epimodel$irm_key_lookup, epimodel$pop_mat, epimodel$index_state_num)

          # get hazards and event rates
          .hazards            <- rep(0, length(.left_endpoints))
          .event_rates        <- rep(0, length(.left_endpoints))

          for(r in 1:length(.hazards)) {

                    .irm <- epimodel$irm[,,.irm_keys[r]]

                    # get the hazard from the appropriate rate matrix
                    .hazards[r] <- .irm[.subj_path[r], .subj_path[r]]

                    # check if the event corresponded to the subject. if so, get the event rate
                    if(epimodel$pop_mat[.right_endpoints[r], "ID"] == subject) {

                              .event_rates[r] <- .irm[.subj_path[r], .subj_path[r + 1]]

                    }
          }

          # calculate subject level likelihood
          subj_likelihood <- log(epimodel$initdist[.subj_path[1]]) + sum(log(.event_rates[.event_rates != 0])) + sum(.hazards * .time_diffs)

          # exponentiate if desired
          if(log == FALSE) {
                    subj_likelihood <- exp(subj_likelihood)
          }

          return(subj_likelihood)
}