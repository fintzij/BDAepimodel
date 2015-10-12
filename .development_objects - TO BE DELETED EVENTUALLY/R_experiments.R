
# Arrays vs. lists --------------------------------------------------------

x <- array(rnorm(1:100000), dim = c(10, 10, 1000))

y <- lapply(seq(dim(x)[3]), function(z) x[,,z])

row.inds <- sample.int(10, 1000, replace = T)
col.inds <- sample.int(10, 1000, replace = T)


mat <- matrix(1:100, nrow = 10)

system.time(replicate(1000, {
          
          # assignment
          x[,,1:dim(x)[3]] <- mat
          
          #retrieval
          x[cbind(row.inds, col.inds, 1:dim(x)[3])]
          
}))

system.time(replicate(1000, {
          
          # assignment
          y[1:dim(x)[3]] <- list(mat)
          
          #retrieval
          z <- array(unlist(y), dim = c(10, 10, length(y)))
          # mapply("[", y, row.inds, col.inds)
          z[cbind(row.inds, col.inds, 1:dim(z)[3])]
          
}))


# Fun with environments ---------------------------------------------------

e <- new.env()
e[["one"]] <- 4
e[["two"]] <- 4
e[["three"]] <- 6
e[["I2"]] <- 7


# Figuring out how to call a function within the environment of another function ----
x <- new.env()
x$a <- 1
x$b <- 1

y <- list(a = 1, b = 2)
z <- list2env(y)

fun1 <- function(.epimodel) {
          .epimodel$config_mat[.epimodel$config_mat[,"ID"] == 1,] <- 999
          .epimodel$log_likelihood <- 361
}

fun2 <- function(epimodel) {
          epimodel$config_mat[epimodel$config_mat[,"ID"] == 1,] <- 999
          return(epimodel$config_mat)
}

fun3 <- function(epimodel) {
          log_likelihood <- 361
          return(log_likelihood)
}

system.time({for(k in 1:10000) fun1(.epimodel)})
system.time({for(k in 1:10000) {epimodel$config_mat <- fun2(epimodel); epimodel$log_likelihood <- 361}})
