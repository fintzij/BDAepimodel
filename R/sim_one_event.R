#' Simulate the next event using Gillespie's next reaction method.
#' 
#' Internal function for \code{simulate_epimodel}, which is its parent 
#' environment. This function uses the rate_mat, t, and epimodel objects in its
#' parent environment, modifying them there.
#' 
#' @param config_list list with three named elements - config, rate_mat, and
#'   dt_mat - where config should correspond to a row in
#'   \code{epimodel$config}.
#' @inheritParams simulate_epimodel
#'   
#' @return updated config_list
#' 
#' @export  
#'

sim_one_event <- function(config_list, epimodel, lump) {
          
          # pass objects to internal objects
          .cur_config <- config_list$config
          .dt_mat <- config_list$dt_mat
          .rate_mat <- config_list$rate_mat
          .t <- .cur_config[,"time"]
          
          # reset dt_mat and rate_mat
          .dt_mat[.dt_mat != Inf]       <- Inf
          .rate_mat[.rate_mat != 0]     <- 0
          
          # get the current configuration and initialize the next configuration
          .subj_config <- .cur_config[,grepl(".X", colnames(.cur_config)), drop = FALSE]
          
          # find the rates for each possible state change
          .rates <- sapply(epimodel$rates, do.call, list(state = .cur_config, params = epimodel$params))
          
          if(lump == TRUE) {
                    for(.r in seq_along(.rates)) {
                              .rates[.r] <- .rates[.r] * sum(.subj_config == which(epimodel$flow[.r, , drop = FALSE] == -1))
                    }
          }
          
          # if not all the rates are 0, sample the next event
          if(!all(.rates == 0)) {
                    
                    if(lump == TRUE) {
                              
                              .rate_mat <- .rates
                              
                              .dt_mat[.rate_mat != 0] <- rexp(n = sum(.rate_mat != 0), rate = .rate_mat[.rate_mat != 0])
                              .event_ID <- which.min(.dt_mat)
                              
                              # get next event - first element corresponds to subject ID, second specifies the transition
                              .next_event <- c(sample.int(n = epimodel$popsize, size = 1, prob = (.subj_config == which(epimodel$flow[.event_ID, , drop = FALSE] == -1))), .event_ID)
                              .next_time <- .t + .dt_mat[.event_ID]
                              
                    } else if(lump == FALSE) {
                              
                              # apply the rates to the appropriate subjects in the rate matrix
                              for(.r in seq_along(.rates)){
                                        
                                        # if the rate is 0, move on to the next iteration
                                        if(.rates[.r] == 0) next
                                        
                                        # otherwise fill in the rates and compute event times
                                        .which_state <- which(epimodel$flow[.r, , drop = FALSE] == -1)
                                        .rate_mat[.subj_config == .which_state, .r] <- .rates[.r]
                                        .dt_mat[.subj_config == .which_state, .r] <- rexp(n = sum(.subj_config == .which_state), rate = .rate_mat[.subj_config == .which_state, .r])
                                        
                              }
                              
                              # get next event - first element corresponds to subject ID, second element specifies the transition
                              .next_event                   <- arrayInd(which.min(.dt_mat), dim(.dt_mat))
                              .next_time                    <- .t + .dt_mat[.next_event[1], .next_event[2]]
                              
                    }
                    
                    # if the next event time is not after tmax, update the vector
                    if(.next_time < max(epimodel$obstimes)){
                              # update the next configuration
                              .cur_config[,"time"]         <- .next_time
                              .cur_config[,"ID"]           <- .next_event[1]
                              .cur_config[,"Event"]        <- .next_event[2]
                              .cur_config[,epimodel$states]<- .cur_config[, epimodel$states] + epimodel$flow[.next_event[2], ]
                              .cur_config[,colnames(.subj_config)[.next_event[1]]] <- which(epimodel$flow[.next_event[2],] == 1)
                              
                              return(list(config = .cur_config, rate_mat = .rate_mat, dt_mat = .dt_mat, keep_going = TRUE))   
                              
                    } else { # return the final observation time
                              .cur_config[,"time"] <- max(epimodel$obstimes)
                              .cur_config[,"ID"] <- 0
                              .cur_config[,"Event"] <- 0
                              
                              return(list(config = .cur_config, rate_mat = .rate_mat, dt_mat = .dt_mat, keep_going = FALSE))
                    }
                    
          } else {
                    .cur_config[,"time"] <- max(epimodel$obstimes)
                    .cur_config[,"ID"] <- 0
                    .cur_config[,"Event"] <- 0
                    
                    return(list(config = .cur_config, rate_mat = .rate_mat, dt_mat = .dt_mat, keep_going = FALSE))
          }

}