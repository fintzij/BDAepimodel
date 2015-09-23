
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
