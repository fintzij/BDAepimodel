#' Simulate the next event using Gillespie's next reaction method.
#'
#' Internal function for \code{simulate_epimodel}, which is its parent
#' environment. This function uses the rate_mat, t, and epimodel objects in its
#' parent environment, modifying them there.
#'
#' @param config_list list with three named elements - config_mat, pop_mat,
#'   rate_mat, and dt_mat.
#' @inheritParams simulate_epimodel
#'
#' @return updated config_list
#'
#' @export
#'

sim_one_event <- function(config_list, epimodel, lump) {

          # pass objects to internal objects
          .pop_config_cur <- config_list$pop_config
          .subj_config_cur <- config_list$subj_config
          .dt_mat <- config_list$dt_mat
          .rate_mat <- config_list$rate_mat
          .t <- .pop_config_cur[,"time"]

          # reset dt_mat and rate_mat
          .dt_mat[.dt_mat != Inf]       <- Inf
          .rate_mat[.rate_mat != 0]     <- 0

          # find the rates for each possible state change
          .rates <- sapply(epimodel$rates, do.call, list(state = .pop_config_cur, params = epimodel$params))

          if(lump == TRUE) {
                    for(.r in seq_along(.rates)) {
                              .rates[.r] <- .rates[.r] * sum(.subj_config_cur == which(epimodel$flow[.r, , drop = FALSE] == -1))
                    }
          }

          # if not all the rates are 0, sample the next event
          if(!all(.rates == 0)) {

                    if(lump == TRUE) {

                              .rate_mat <- .rates

                              .dt_mat[.rate_mat != 0] <- rexp(n = sum(.rate_mat != 0), rate = .rate_mat[.rate_mat != 0])
                              .event_ID <- which.min(.dt_mat)

                              # get next event - first element corresponds to subject ID, second specifies the transition
                              .next_event <- c(sample.int(n = epimodel$popsize, size = 1, prob = (.subj_config_cur == which(epimodel$flow[.event_ID, , drop = FALSE] == -1))), .event_ID)
                              .next_time <- .t + .dt_mat[.event_ID]

                    } else if(lump == FALSE) {

                              # apply the rates to the appropriate subjects in the rate matrix
                              for(.r in seq_along(.rates)){

                                        # if the rate is 0, move on to the next iteration
                                        if(.rates[.r] == 0) next

                                        # otherwise fill in the rates and compute event times
                                        .which_state <- which(epimodel$flow[.r, , drop = FALSE] == -1)
                                        .rate_mat[.subj_config_cur == .which_state, .r] <- .rates[.r]
                                        .dt_mat[.subj_config_cur == .which_state, .r] <- rexp(n = sum(.subj_config_cur == .which_state), rate = .rate_mat[.subj_config_cur == .which_state, .r])

                              }

                              # get next event - first element corresponds to subject ID, second element specifies the transition
                              .next_event                   <- arrayInd(which.min(.dt_mat), dim(.dt_mat))
                              .next_time                    <- .t + .dt_mat[.next_event[1], .next_event[2]]

                    }

                    # if the next event time is not after tmax, update the vector
                    if(.next_time < max(epimodel$obstimes)){
                              # update the next configuration
                              .pop_config_cur[,"time"]         <- .next_time
                              .pop_config_cur[,"ID"]           <- .next_event[1]
                              .pop_config_cur[,"Event"]        <- .next_event[2]
                              .pop_config_cur[,epimodel$states]<- .pop_config_cur[, epimodel$states] + epimodel$flow[.next_event[2], ]
                              .subj_config_cur[,colnames(.subj_config_cur)[.next_event[1]]] <- which(epimodel$flow[.next_event[2],] == 1)

                              return(list(pop_config = .pop_config_cur, subj_config = .subj_config_cur, rate_mat = .rate_mat, dt_mat = .dt_mat, keep_going = TRUE))

                    } else { # return the final observation time
                              .pop_config_cur[,"time"] <- max(epimodel$obstimes)
                              .pop_config_cur[,"ID"] <- 0
                              .pop_config_cur[,"Event"] <- 0

                              return(list(pop_config = .pop_config_cur, subj_config = .subj_config_cur, rate_mat = .rate_mat, dt_mat = .dt_mat, keep_going = FALSE))
                    }

          } else {
                    .pop_config_cur[,"time"] <- max(epimodel$obstimes)
                    .pop_config_cur[,"ID"] <- 0
                    .pop_config_cur[,"Event"] <- 0

                    return(list(pop_config = .pop_config_cur, subj_config = .subj_config_cur, rate_mat = .rate_mat, dt_mat = .dt_mat, keep_going = FALSE))
          }
}