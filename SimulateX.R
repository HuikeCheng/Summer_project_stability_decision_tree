transformX <- function(x, association, sd = 0.1, width = 1) {
  Functions <- list(Exponential = function(x){exp(x)},
                    Quadratic = function(x){x*x},
                    Log = function(x){log(x)},
                    Sqrt = function(x){sqrt(x)},
                    Sigmoidal = function(x){1/(1+exp(-2*x))},
                    Sine = function(x){2*sin(x)},
                    Heteroscedastic = function(x){x*runif(n = length(x), min=-3, max=3)*rnorm(n = length(x), sd = sd)},
                    Absolute = function(x){abs(x)},
                    Sawtooth = function(x){gsignal::sawtooth(x, width = width)},
                    SquareWave = function(x){gsignal::square(x)},
                    Cubic = function(x) {x^3},
                    Indicator = function(x) {as.numeric(x>0)},
                    Cosine = function(x) {2*cos(x)})
  if (association == "Linear") {
    x <- x
  } else if (!(association %in% names(Functions))) {
    stop("Not recognised")
  } else {
    tmp <- which(names(Functions) == association)
    x <- Functions[[tmp]](x)
  }
  if (is.matrix(x) == FALSE) {
    x <- as.matrix(x)
  }
  return(x)
}

# combineX <- function(X, type) {
#   if (type == "Addition") {
#     beta <- rep(1, ncol(X))
#     out <- X%*%beta
#   } else if (type == "Multiplication") {
#     out <- X[,1]
#     for (i in 2:ncol(X)) {
#       out <- out*X[,i]
#     }
#   }
#   return(out)
# }

HugeAdjacency <- function(pk = 10, topology = "random", nu = 0.1, ...) {
  # Storing extra arguments
  extra_args <- list(...)
  
  # Extracting relevant extra arguments
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = huge::huge.generator)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("n", "d", "prob", "graph", "verbose")]
  
  # Running simulation model
  mymodel <- do.call(huge::huge.generator, args = c(
    list(
      n = 2, d = sum(pk), prob = nu,
      graph = topology, verbose = FALSE
    ),
    tmp_extra_args
  ))
  theta <- as.matrix(mymodel$theta)
  
  # Re-organising the variables to avoid having centrality related to variable ID (e.g. for scale-free models)
  ids <- sample(ncol(theta))
  theta <- theta[ids, ids]
  
  return(theta)
}


# library(ggplot2)
# library(ggpubr)
# association <- c("Linear","Cubic", "Sigmoidal", "Exponential", "Sawtooth", "SquareWave", 
#                  "Indicator", "Sine", "Cosine","Quadratic", "Absolute", "Heteroscedastic")
# out <- vector("list", length = length(association))
# x <- rnorm(1000)
# for (i in seq_along(association)) {
#   y <- SimulateAssoc(x, association=association[i], sd=0.1, width=1)
#   out[[i]] <- ggplot(data.frame(x=x, y=y), aes(x,y)) +
#     geom_point() +
#     ggtitle(association[i])
# }
# 
# ggarrange(out[[1]], out[[2]], out[[3]], out[[4]], out[[5]], out[[6]], out[[7]],
#           out[[8]], out[[9]], out[[10]], out[[11]], out[[12]])


# x <- rnorm(1000)
# y <- SimulateAssoc1(x,"Indicator")
# ggplot(data.frame(x=x, y=y), aes(x,y)) +
#      geom_point() +
#      ggtitle("Indicator function")

# n <- 1000 # number of observations 
# pk <- 100 #number of predictors 
# 
# X <- matrix(rnorm(n*pk),n,pk)
# 
# R <- 1:3
# S <- 4:6 
# 
# alpha <- 10*runif(3)
# beta <- 10*runif(3) 
# 
# mu.1 <- X[,R]%*%alpha #X1+X2+X3
# 
# mu.2 <- X[,S]%*%beta #X4+X5+X6
# 
# Y <- 5+ 3*(mu.1>0)*5*(mu.2>0) + rnorm(n,0,1) 
# 
# plot3d(mu.1,mu.2,Y)