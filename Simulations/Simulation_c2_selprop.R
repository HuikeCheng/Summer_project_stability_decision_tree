######## pass in arguments
args=commandArgs(trailingOnly=TRUE)
#args <- c("100", "1000", "500", "0.01", "5")
numrep <- as.numeric(args[1])
n <- as.numeric(args[2])
pk <- as.numeric(args[3])
ev_xy <- as.numeric(args[4])
nchunks <- as.numeric(args[5])

######## create folder
folder <- paste(paste0("SelpropOutputs/", "Complex2"), n, pk, ev_xy, sep = "_")
if (file.exists(folder) == FALSE) {
  dir.create(folder)
}
print(folder)

######### load packages
library(rpart)
library(sharp)
#library(MASS)
library(rlist)
library(tidyverse)
library(parallel)
#library(pbapply)
library(ggplot2)
source("../cart1.R")
source("../cart2.R")
source("../Functions.R")

####################### Set up simulation study ################################
Simulation_study <- function(seed, n, pk, ev_xy, multivariate_normal = TRUE) {
  ##### set seed
  set.seed(seed)
  ##### simulate data
  # simulate xdata
  if (multivariate_normal == TRUE) {
    xsimul <- SimulateGraphical(
      n = n, pk = pk, theta = NULL,
      implementation = HugeAdjacency, topology = "random",
      nu_within = 0, nu_between = 0, nu_mat = NULL,
      v_within = 0, v_between = 0,
      v_sign = c(-1, 1), continuous = TRUE,
      pd_strategy = "diagonally_dominant", ev_xx = NULL, scale_ev = TRUE,
      u_list = c(1e-10, 1), tol = .Machine$double.eps^0.25,
      scale = TRUE, output_matrices = FALSE
    )
    xdata <- xsimul$data
    print("mn")
  } else {
    xdata <- matrix(rnorm(n*pk),n,pk)
    colnames(xdata) <- paste0("var", 1:pk)
  }
  # simulate ydata
  ypred <- 0.25*exp(4*xdata[,1]) + 4/(1+exp(-20*(xdata[,2]-0.5))) + xdata[,3]
  
  sigma <- sqrt((1 - ev_xy) / ev_xy * stats::var(ypred))
  ydata <- stats::rnorm(n = n, mean = ypred, sd = sigma)
  ydata <- as.matrix(ydata)
  
  ######### run stab1
  scart1 <- sharp::VariableSelection(
    xdata = xdata, ydata = ydata,
    implementation = CART1,
    family = "gaussian",
    maxsurrogate = 0)
  
  selprops <- scart1$selprop[ArgmaxId(scart1)[1],]
  selthr <- Argmax(scart1)[2]
  
  ####### output
  return(list(selprops = selprops, selthr = selthr))
}

# start.time <- Sys.time()
# a <- Simulation_study(seed = 1, n = 500, pk = 1000, ev_xy = 0.2)
# end.time <- Sys.time()
# end.time - start.time

########### Parallelise
no_cores <- nchunks
cl <- makeCluster(no_cores)

clusterEvalQ(cl, library(rpart))
clusterEvalQ(cl, library(sharp))
clusterEvalQ(cl, library(fake))
clusterEvalQ(cl, library(rlist))
clusterExport(cl, c("numrep", "n", "pk", "ev_xy", "HugeAdjacency",
                    "CART1", "CART2", "getMetrics", "Simulation_study"))

out <- parLapply(1:numrep,
                 function(i) {Simulation_study(seed = i, n = n, pk = pk, ev_xy = ev_xy)},
                 cl = cl)

stopCluster(cl)

########### cleanup output
Selprops <- sapply(out, "[", 1)
Selprops <- list.cbind(Selprops)

Selthrs <- sapply(out, "[", 2)
Selthrs <- list.cbind(Selthrs)

theta <- rep(0, pk)
names(theta) <- paste0("var", 1:pk)
theta[1:3] <- 1

########### save output
saveRDS(Selprops, file = paste(folder,"Selprops.rds", sep = "/"))
saveRDS(Selthrs, file = paste(folder,"Selthrs.rds", sep = "/"))
saveRDS(theta, file = paste(folder,"Theta.rds", sep = "/"))