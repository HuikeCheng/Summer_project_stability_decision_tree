SimulateAssoc <- function(x, association, sd = 0.1, width = 1) {
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
<<<<<<< HEAD
  } else if (association == "Exponential") {
    x <- exp(x)
  } else if (association == "Log") {
    x <- log(x)
  } else if (association == "Quadratic") {
    x <- x*x
  } else if (association == "Sqrt") {
    x <- sqrt(x)
  } else if (association == "Sigmoidal") {
    x <- 1/(1+exp(-2*x))
  } else if (association == "Sine") {
    x <- 2*sin(x)
  } else if (association == "Heteroscedastic") {
    x <- x*runif(n = length(x), min=-3, max=3)*rnorm(n = length(x), sd = sd)
  } else if (association == "Absolute") {
    x <- abs(x)
  } else if (association == "Threshold") {
    # Simple threshold effect
    if (!is.null(threshold)) {
      x[x<=threshold] <- min(x[x>=threshold])
    }
=======
  } else if (!(association %in% names(Functions))) {
    stop("Not recognised")
>>>>>>> c8235c1f41c8b8394c949083a96e6e51b66c3958
  } else {
    tmp <- which(names(Functions) == association)
    x <- Functions[[tmp]](x)
  }
  if (is.matrix(x) == FALSE) {
    x <- as.matrix(x)
  }
  return(x)
}

CombineX <- function(x, type) {
  if (type == "Addition") {
    beta <- rep(1, ncol(x))
    out <- x%*%beta
  } else if (type == "Multiplication") {
    out <- x[,1]
    for (i in 2:ncol(x)) {
      out <- out*x[,i]
    }
  }
  return(out)
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