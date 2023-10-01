############
args <- c("10", "5", "500", "1000", "0.2", "5")
numrep <- as.numeric(args[1])
height <- as.numeric(args[2])
n <- as.numeric(args[3])
pk <- as.numeric(args[4])
ev_xy <- as.numeric(args[5])
nchunks <- as.numeric(args[6])

######### load packages
library(rpart)
library(sharp)
library(fake)
library(rlist)
library(tidyverse)
library(data.tree)
library(parallel)
library(pbapply)
library(ggplot2)
source("cart1.R")
source("cart2.R")
source("Functions.R")
source("SimulateTree2.R")

##############################
Computation_time <- function(seed, height, n, pk, ev_xy) {
  # set seed
  set.seed(seed)
  # generate data
  simul <- SimulateTree(height = height, n = n, pk = pk, ev_xy = ev_xy)
  
  # set up output
  mytime <- vector("list", length = 6)
  names(mytime) <- c("CART", "sCART1", "sCART1t","sCART2", "sCART2t","sLASSO")
  
  ######## run CART
  mydata <- data.frame(ydata = simul$ydata[, 1], simul$xdata)
  
  start.time <- Sys.time()
  cart <- rpart(ydata ~ ., mydata, method = "anova", maxsurrogate = 0)
  end.time <- Sys.time()
  
  time.taken <- difftime(end.time, start.time, units = "secs")
  mytime$CART <- time.taken
  
  ######### run stab1
  start.time <- Sys.time()
  scart1 <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = CART1,
    family = "gaussian",
    maxsurrogate = 0)
  end.time <- Sys.time()
  
  time.taken <- difftime(end.time, start.time, units = "secs")
  mytime$sCART1 <- time.taken
  
  ######### run stab1_k=1000
  start.time <- Sys.time()
  scart1 <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = CART1,
    family = "gaussian",
    maxsurrogate = 0,
    K=500)
  end.time <- Sys.time()
  
  time.taken <- difftime(end.time, start.time, units = "secs")
  mytime$sCART1t <- time.taken
  
  ####### run stab2
  Lambda <- cart$cptable[,1]
  
  start.time <- Sys.time()
  scart2 <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = CART2,
    family = "gaussian",
    Lambda = Lambda,
    maxsurrogate = 0)
  end.time <- Sys.time()
  
  time.taken <- difftime(end.time, start.time, units = "secs")
  mytime$sCART2 <- time.taken
  
  ####### run stab2, k=1000
  start.time <- Sys.time()
  scart2 <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = CART2,
    family = "gaussian",
    Lambda = Lambda,
    maxsurrogate = 0,
    K=500)
  end.time <- Sys.time()
  
  time.taken <- difftime(end.time, start.time, units = "secs")
  mytime$sCART2t <- time.taken
  
  ####### run Stability-lasso
  start.time <- Sys.time()
  slasso <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = PenalisedRegression,
    family = "gaussian")
  end.time <- Sys.time()
  
  time.taken <- difftime(end.time, start.time, units = "secs")
  mytime$sLASSO <- time.taken
  
  ####### Cleanup output
  mytime <- as.data.frame(list.rbind(mytime))
  mytime$model <- rownames(mytime)
  
  ####### output
  return(mytime)
}

########### Parallelise
numrep <- 10
no_cores <- nchunks
cl <- makeCluster(no_cores)

clusterEvalQ(cl, library(rpart))
clusterEvalQ(cl, library(sharp))
clusterEvalQ(cl, library(fake))
clusterEvalQ(cl, library(rlist))
clusterEvalQ(cl, library(data.tree))
clusterExport(cl, c("numrep", "height", "n", "pk", "ev_xy", "HugeAdjacency",
                    "reorderCV", "CART1", "CART2", "Computation_time",
                    "SimulateTree", "Bin2Dec"))

out <- pblapply(1:numrep,
                function(i) {Computation_time(seed = i, height = height, 
                                              n = n, pk = pk, ev_xy = ev_xy)},
                cl = cl)

stopCluster(cl)

############ cleanup output
#mytime2 <- list.rbind(out)
#mytime2$dim <- rep("n = 1000, p = 500", 60)

#mytime <- list.rbind(out)
#mytime$dim <- rep("n = 500, p = 1000", 60)

#df <- rbind(mytime, mytime2)

#saveRDS(df, "Kdiff.rds")
