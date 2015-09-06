# #' Initializes a matrix of subject level trajectories.
# #' 
# #' @param epimodel bookkeeping object
# #' @inheritParams simulate_epimodel
# #'   
# #' @return matrix of subject level trajectories, where each subject has two
# #'   rows, one for time zero and one for tmax.
# #'   
# init_subj_mat <- function(epimodel, init_state){
#          
#           popsize = sum(init_state)
#           
#           # set up matrix
#           subj_mat <- matrix(0, nrow = 2*epimodel$popsize, ncol = 2 + length(epimodel$states))
#           colnames(subj_mat) <- c("ID", "time", epimodel$states)
#           
#           # individual IDs and times - n.b. tmax = max(obstimes)
#           subj_mat[, "ID"] <- rep(1:popsize, each = 2)
#           subj_mat[, "time"] <- rep(c(0, max(epimodel$obstimes)), epimodel$popsize)
#           
#           # set initial and final statuses according to init_state
#           IDs <- 1:popsize
#           for(t in seq_along(init_state))  {
#                     if(init_state[t] == 0) next
#                     selected <- IDs[1:init_state[t]]
#                     inds <- which(subj_mat[ ,"ID"] %in% selected)[seq(1, sum(subj_mat[ ,"ID"] %in% selected), by = 2)]
#                     subj_mat[inds, names(init_state)[t]] <- 1
#                     IDs <- IDs[-c(1:init_state[t])]
#           }
#           
#           return(subj_mat)
# }