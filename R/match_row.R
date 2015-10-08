#' Returns either a logical vector of indices for which row(s) in a table match 
#' a given vector or rows in a matrix.
#' 
#' @param x vector or matrix
#' @param table table in which to search for matches
#' @param return_value one of "index", "logical", or "vector" to return the index
#'   of the matching row (possibly of length 0), a logical value indicating
#'   whether there is a matching row, or a vector of logicals indicating a match
#'   or not.
#'   
#' @return vector of logicals or indices
#' @export
#' 
match_row <- function(x, table, return_value) {
          
          if(length(x) == 1L) {
                    .inds <- .Internal(which(table[, 1, drop = FALSE] == as.numeric(x)))
                    
          } else {
                    # initialize vector of TRUE values
                    .inds <- 1:nrow(table)
                    
                    # find matches
                    for(j in 1:length(x)) {
                              .inds <- .inds[table[.inds, j, drop = FALSE] == x[j]]
                    }
                    
          }

          
          if(return_value == "logical") {
                    return(any(.inds)) 
                    
          } else if(return_value == "index") {
                    return(.inds)
                    
          } else if(return_value == "vector") {
                    return(replace(logical(length = nrow(table)), .inds, T))
                    
          } 
}
