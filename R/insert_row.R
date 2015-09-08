#' Generic function to insert a row into a matrix.
#'
#' @param mat matrix into which row should be inserted
#' @param row_num row number for insertion; row currently at that number is pushed down.
#' @param vec vector to be inserted
#'
#' @return matrix with the inserted row. 
#'
#' @examples x <- matrix(1:9, nrow = 3)
#' insert_row(x, 2, c(5, 27, 87))
#' 
insert_row <- function(mat, row_num, vec) {
          
          nrows <- nrow(mat)
          ncols <- ncol(mat)

          if (row_num == 1) {
                    new_mat <- rbind(matrix(vec, ncol = ncols), mat, deparse.level = 0)
          }
          else if (row_num == nrows + 1) {
                    new_mat <- rbind(mat, matrix(vec, ncol = ncols), deparse.level = 0)
          }
          else {
                    new_mat <- rbind(mat[1:(row_num - 1), , drop = FALSE],
                                 matrix(vec, ncol = ncols), mat[row_num:nrows, , drop = FALSE])
          }
          return(new_mat)
}