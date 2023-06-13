SimulateAssoc <- function(x, association, sd = 0.1, threshold = NULL) {
  if (association == "Linear") {
    x <- x
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
  }
  
  # Simple threshold effect
  if (!is.null(threshold)) {
    x[x<=threshold] <- min(x[x>=threshold])
  }
  return(x)
}


# library(ggplot2)
# library(ggpubr)
# association <- "Linear"
# set.seed(1)
# x <- rnorm(1000)
# x <- abs(x)
# y <- SimulateAssoc(x, association = association, threshold = -1)
# p8 <- ggplot(data.frame(x=x, y=y), aes(x, y)) + 
#   geom_point() +
#   ggtitle("Linear with simple threshold effect")
# 
# ggarrange(p1,p2,p3,p4,p5,p6,p7,p8)


#SimulateInteraction
#number of pairs
#invThreshold
#vThreshold

# x <- 1
# a <- "2*sin(x)"
# eval(parse(a))
# 
# b <- function(x, assoc){
#   y = assoc(x)
#   return(y)
#   }
# 
# b(1, assoc=function(x) {2*sin(x)})

# z <- runif(10, -0.5, 0.5) + 4*rbinom(10,1,0.5)-2
# w <- runif(10, -0.5, 0.5) + 4*rbinom(10,1,0.5)-2
# x <- cos(pi/128)*z+sin(pi/128)*w
