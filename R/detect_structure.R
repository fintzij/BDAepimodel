#' Detects whether there is structure in the stochastic epidemic model that can 
#' be leveraged.
#' 
#' Detects whether the model has absorbing states (exit from such states is 
#' impossible), and whether the model is progressive. If there is an absorbing 
#' state, no additional path sampling is needed once a subject has entered that 
#' state. If the model is progressive, it only moves forward through states and 
#' never returns to a previous state. Therefore, if the model is progressive and
#' has an absorbing state, only the transition times in intervals where the
#' state has changed need to be resampled. 
#' 
#' @param epimodel
#'   
#' @return objects in the epimodel environment: absorbing states (or FALSE if 
#'   none), and with TRUE/FALSE for whether the model is progressive.
#' @export  
#'   

detect_structure <- function(epimodel) {
          
          # detect whether each state is absorbing 
          absorbing_states <- rep(FALSE, epimodel$num_states)
          
          for(j in 1:epimodel$num_states) {
                    absorbing_states[j] <-  ifelse((1 %in% epimodel$flow[,j]) & (!-1 %in% epimodel$flow[,j]), TRUE, FALSE)
          }
          
          epimodel$absorbing_states <- ifelse(any(absorbing_states), which(absorbing_states), FALSE)
          
          # detect whether the model is progressive 
          dag <- matrix(0, nrow = epimodel$num_states, ncol = epimodel$num_states)
          
          for(j in 1:nrow(epimodel$flow)) {
                    dag[which(epimodel$flow[j,, drop = FALSE] == -1), which(epimodel$flow[j,, drop = FALSE] == 1)] <- 1
          }
          
          epimodel$progressive <- ifelse(any(dag[lower.tri(dag)] != 0), FALSE, TRUE)
          
}