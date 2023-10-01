######## pass in arguments
args=commandArgs(trailingOnly=TRUE)
#args <- c("100", "Linear", "1000", "500", "0.2", "2", "2", "50")
numrep <- as.numeric(args[1])
association <- as.character(args[2])
n <- as.numeric(args[3])
pk <- as.numeric(args[4])
ev_xy <- as.numeric(args[5])
numpair <- as.numeric(args[6])
numvar <- as.numeric(args[7])
nchunks <- as.numeric(args[8])

######## create folder
folder <- paste(paste0("ReMult4Outputs/", "Mult"), association, n, pk, ev_xy, numpair, numvar, sep = "_")
if (file.exists(folder) == FALSE) {
  dir.create(folder)
}
print(folder)
#####################
library(rpart)
library(sharp)
library(fake)
library(rlist)
library(tidyverse)
library(parallel)
#library(pbapply)
library(ggplot2)
source("../Functions.R")
source("../cart1.R")
source("../cart2.R")
source("../SimulateX.R")

#####################
Simulation_study <- function(seed, n, pk, association, ev_xy, numpair, numvar, multivariate_normal = TRUE) {
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
  if (numpair == 1) {
    if (numvar == 2) {
      mu <- xdata[,1]*xdata[,2]
    } else if (numvar > 2) {
      mu <- xdata[,1]
      for (i in 2:numvar) {
        mu <- mu*xdata[,i]
      }
    }
    ypred <- transformX(mu, association)
  } else if (numpair == 2) {
    mu1 <- xdata[,1]*xdata[,2]
    mu2 <- xdata[,3]*xdata[,4]
    ypred <- transformX(mu1, association) * transformX(mu2, association)
  }
  
  sigma <- sqrt((1 - ev_xy) / ev_xy * stats::var(ypred))
  ydata <- stats::rnorm(n = n, mean = ypred, sd = sigma)
  ydata <- as.matrix(ydata)
  
  # theta
  theta <- rep(0, pk)
  names(theta) <- colnames(xdata)
  theta[1:(numpair*numvar)] <- 1
  
  # set up output
  metrics <- vector("list", length = 4)
  sel <- vector("list", length = 4)
  names(metrics) <- c("CART", "sCART1", "sCART2", "sLASSO")
  names(sel) <- c("CART", "sCART1", "sCART2", "sLASSO")
  
  ######## run CART
  mydata <- data.frame(ydata = ydata[, 1], xdata)
  cart <- rpart(ydata ~ ., mydata, method = "anova", maxsurrogate = 0)
  
  myvars <- unique(rownames(cart$splits))
  selvars <- rep(0, pk)
  names(selvars) <- colnames(xdata)
  selvars[which(names(selvars) %in% myvars)] <- 1
  metrics$CART <- getMetrics(selvars = selvars, theta = theta)
  sel$CART <- selvars
  
  ######### run stab1
  scart1 <- sharp::VariableSelection(
    xdata = xdata, ydata = ydata,
    implementation = CART1,
    family = "gaussian",
    maxsurrogate = 0)
  
  selvars <- SelectedVariables(scart1)
  metrics$sCART1 <- getMetrics(selvars = selvars, theta = theta)
  sel$sCART1 <- selvars
  
  ####### run stab2
  Lambda <- cart$cptable[,1]
  scart2 <- sharp::VariableSelection(
    xdata = xdata, ydata = ydata,
    implementation = CART2,
    family = "gaussian",
    Lambda = Lambda,
    maxsurrogate = 0)
  
  selvars <- SelectedVariables(scart2)
  metrics$sCART2 <- getMetrics(selvars = selvars, theta = theta)
  sel$sCART2 <- selvars
  
  ####### run Stability-lasso
  slasso <- sharp::VariableSelection(
    xdata = xdata, ydata = ydata,
    implementation = PenalisedRegression,
    family = "gaussian")
  
  selvars <- SelectedVariables(slasso)
  metrics$sLASSO <- getMetrics(selvars = selvars, theta = theta)
  sel$sLASSO <- selvars
  
  ####### Cleanup output
  metrics <- as.data.frame(list.rbind(metrics))
  metrics$model <- rownames(metrics)
  
  sel <- as.data.frame(list.rbind(sel))
  sel$model <- rownames(sel)
  
  ####### output
  return(list(metrics = metrics, sel = sel))
}

#a <- Simulation_study(seed = 1, n = 100, pk = 500, association = "Linear", ev_xy = 0.2, numpair = 2, numvar = 2)

#########################
no_cores <- nchunks
cl <- makeCluster(no_cores)

clusterEvalQ(cl, library(rpart))
clusterEvalQ(cl, library(sharp))
clusterEvalQ(cl, library(fake))
clusterEvalQ(cl, library(rlist))
clusterExport(cl, c("n", "pk", "association","ev_xy", "numpair", "numvar", "numrep",
                    "CART1", "CART2", "getMetrics", "Simulation_study",
                    "HugeAdjacency","transformX"))

out <- parLapply(1:numrep,
                function(i) {Simulation_study(seed = i, n = n, pk = pk, association = association, 
                                              ev_xy = ev_xy, numpair = numpair, numvar = numvar)},
                cl = cl)

stopCluster(cl)

########### cleanup output
Metrics <- sapply(out, "[", 1)
Metrics <- list.rbind(Metrics)
Metrics$model <- factor(Metrics$model, levels = unique(Metrics$model))

Sels <- sapply(out, "[", 2)
Sels <- list.rbind(Sels)

theta <- rep(0, pk)
names(theta) <- paste0("var", 1:pk)
theta[1:(numpair*numvar)] <- 1

########### save output
saveRDS(Metrics, file = paste(folder,"Metrics.rds", sep = "/"))
saveRDS(Sels, file = paste(folder,"Sels.rds", sep = "/"))
saveRDS(theta, file = paste(folder,"Theta.rds", sep = "/"))
