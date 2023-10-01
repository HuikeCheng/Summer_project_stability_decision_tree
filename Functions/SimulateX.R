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