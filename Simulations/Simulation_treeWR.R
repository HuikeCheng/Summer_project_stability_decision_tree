######## pass in arguments
args=commandArgs(trailingOnly=TRUE)
#args <- c("10", "5", "500", "50", "0.1", "5")
numrep <- as.numeric(args[1])
height <- as.numeric(args[2])
n <- as.numeric(args[3])
pk <- as.numeric(args[4])
ev_xy <- as.numeric(args[5])
nchunks <- as.numeric(args[6])

######## create folder
folder <- paste(paste0("OutputsV3/TreeWR"), height, n, pk, ev_xy, sep = "_")
if (file.exists(folder) == FALSE) {
  dir.create(folder)
}
print(folder)
#####################
library(rpart)
library(partykit)
library(sharp)
library(fake)
library(rlist)
library(tidyverse)
library(randomForest)
library(Boruta)
library(parallel)
library(ggplot2)
library(xgboost)
library(rfvimptest)
library(data.tree)
source("Functions/Functions.R")
source("Functions/cart1.R")
source("Functions/cart2.R")
source("Functions/cart3.R")
source("Functions/SimulateX.R")

#####################
time1 <- Sys.time()

#####################
Simulation_study <- function(seed, height, n, pk, ev_xy) {
  # set seed
  set.seed(seed)
  # generate data
  simul <- SimulateTree(height = height, n = n, pk = pk, ev_xy = ev_xy, repeated_pred = TRUE)
  ydata <- simul$ydata
  xdata <- simul$xdata
  # theta
  theta <- simul$theta
  
  # check correlation of xdata
  xCors <- abs(cor(xdata)[upper.tri(cor(xdata), diag = FALSE)])
  if(any(xCors > 0.5)) {
    print("There are highly correlated variables")
  }
  
  # set up output
  metrics <- vector("list", length = 11)
  sel <- vector("list", length = 11)
  names(metrics) <- c("CART", "sCART1", "sCART2", "sCART3", "sLASSO", "ciTree", "RF",
                      "RF_B_CT", "RF_B_C", "XGB_B_CT", "XGB_B_C")
  names(sel) <- c("CART", "sCART1", "sCART2", "sCART3", "sLASSO", "ciTree", "RF",
                  "RF_B_CT", "RF_B_C", "XGB_B_CT", "XGB_B_C")
  
  selprops <- vector("list", length = 4)
  selthr <- vector("list", length = 4)
  names(selprops) <- c("sCART1", "sCART2", "sCART3", "sLASSO")
  names(selthr) <- c("sCART1", "sCART2", "sCART3", "sLASSO")
  
  ######## run CART
  mydata <- data.frame(ydata = ydata[, 1], xdata)
  cart <- rpart(ydata ~ ., mydata, method = "anova", maxsurrogate = 0)
  
  myvars <- unique(rownames(cart$splits))
  selvars <- rep(0, pk)
  names(selvars) <- colnames(xdata)
  selvars[which(names(selvars) %in% myvars)] <- 1
  metrics$CART <- getMetrics(selvars = selvars, theta = theta)
  sel$CART <- selvars
  
  ######### run scart1
  scart1 <- sharp::VariableSelection(
    xdata = xdata, ydata = ydata,
    implementation = CART1,
    family = "gaussian",
    Lambda = cbind(seq(1, min(nrow(xdata) / 2, 100))),
    maxsurrogate = 0)
  
  selvars <- SelectedVariables(scart1)
  names(selvars) <- names(theta)
  metrics$sCART1 <- getMetrics(selvars = selvars, theta = theta)
  sel$sCART1 <- selvars
  
  selprops$sCART1 <- scart1$selprop[ArgmaxId(scart1)[1],]
  selthr$sCART1 <- Argmax(scart1)[2]
  
  ####### run scart2
  Lambda <- cart$cptable[,1]
  scart2 <- sharp::VariableSelection(
    xdata = xdata, ydata = ydata,
    implementation = CART2,
    family = "gaussian",
    Lambda = Lambda,
    maxsurrogate = 0)
  
  selvars <- SelectedVariables(scart2)
  names(selvars) <- names(theta)
  metrics$sCART2 <- getMetrics(selvars = selvars, theta = theta)
  sel$sCART2 <- selvars
  
  selprops$sCART2 <- scart2$selprop[ArgmaxId(scart2)[1],]
  selthr$sCART2 <- Argmax(scart2)[2]
  
  ####### run scart3
  scart3 <- sharp::VariableSelection(
    xdata = xdata, ydata = ydata,
    Lambda = cbind(seq(0, 0.2, 0.005)),
    implementation = CART3)
  
  selvars <- SelectedVariables(scart3)
  names(selvars) <- names(theta)
  metrics$sCART3 <- getMetrics(selvars = selvars, theta = theta)
  sel$sCART3 <- selvars
  
  selprops$sCART3 <- scart3$selprop[ArgmaxId(scart3)[1],]
  selthr$sCART3 <- Argmax(scart3)[2]
  
  ####### run Stability-lasso
  slasso <- sharp::VariableSelection(
    xdata = xdata, ydata = ydata,
    implementation = PenalisedRegression,
    family = "gaussian")
  
  selvars <- SelectedVariables(slasso)
  metrics$sLASSO <- getMetrics(selvars = selvars, theta = theta)
  sel$sLASSO <- selvars
  
  selprops$sLASSO <- slasso$selprop[ArgmaxId(slasso)[1],]
  selthr$sLASSO <- Argmax(slasso)[2]
  
  ####### run ctree
  citree <- partykit::ctree(ydata~., data = mydata)
  pvals <- nodeapply(citree, ids = nodeids(citree), FUN = function(n) info_node(n)$p.value)
  if (list(NULL) %in% pvals) {
    pvals <- Filter(Negate(is.null), pvals)
  }
  
  if (any(pvals < 0.05)) {
    selvars <- pvals[which(pvals < 0.05)]
    myvars <- sapply(selvars, function(x){names(x[1])})
  }
  
  selvars <- rep(0, pk)
  names(selvars) <- colnames(xdata)
  selvars[which(names(selvars) %in% myvars)] <- 1
  metrics$ciTree <- getMetrics(selvars = selvars, theta = theta)
  sel$ciTree <- selvars
  
  ####### run RF
  rf <- randomForest(x = xdata, y = ydata)
  rf_importance <- rf$importance[order(rf$importance, decreasing = TRUE),]
  myvars <- names(rf_importance)[1:numvar]
  selvars <- rep(0, pk)
  names(selvars) <- colnames(xdata)
  selvars[which(names(selvars) %in% myvars)] <- 1
  metrics$RF <- getMetrics(selvars = selvars, theta = theta)
  sel$RF <- selvars
  
  ####### run RF_B
  rf_b <- Boruta(ydata~.,data=mydata,doTrace=2, getImp = getImpLegacyRfZ)
  myvars <- names(rf_b$finalDecision)[which(rf_b$finalDecision != "Rejected")]
  selvars <- rep(0, pk)
  names(selvars) <- colnames(xdata)
  selvars[which(names(selvars) %in% myvars)] <- 1
  metrics$RF_B_CT <- getMetrics(selvars = selvars, theta = theta)
  sel$RF_B_CT <- selvars
  
  myvars <- names(rf_b$finalDecision)[which(rf_b$finalDecision == "Confirmed")]
  selvars <- rep(0, pk)
  names(selvars) <- colnames(xdata)
  selvars[which(names(selvars) %in% myvars)] <- 1
  metrics$RF_B_C <- getMetrics(selvars = selvars, theta = theta)
  sel$RF_B_C <- selvars
  
  ####### XGB_B
  xgb_b <- Boruta(ydata~.,data=mydata,doTrace=2, getImp = getImpXgboost, nthread=1)
  myvars <- names(xgb_b$finalDecision)[which(xgb_b$finalDecision != "Rejected")]
  selvars <- rep(0, pk)
  names(selvars) <- colnames(xdata)
  selvars[which(names(selvars) %in% myvars)] <- 1
  metrics$XGB_B_CT <- getMetrics(selvars = selvars, theta = theta)
  sel$XGB_B_CT <- selvars
  
  myvars <- names(xgb_b$finalDecision)[which(xgb_b$finalDecision == "Confirmed")]
  selvars <- rep(0, pk)
  names(selvars) <- colnames(xdata)
  selvars[which(names(selvars) %in% myvars)] <- 1
  metrics$XGB_B_C <- getMetrics(selvars = selvars, theta = theta)
  sel$XGB_B_C <- selvars
  
  ####### Cleanup output
  metrics <- as.data.frame(list.rbind(metrics))
  metrics$model <- rownames(metrics)
  
  sel <- as.data.frame(list.rbind(sel))
  sel$model <- rownames(sel)
  
  selprops <- as.data.frame(list.rbind(selprops))
  selprops$model <- rownames(selprops)
  
  selthr <- as.data.frame(list.rbind(selthr))
  selthr$model <- rownames(selthr)
  colnames(selthr)[1] <- "selthr"
  
  ####### output
  return(list(metrics = metrics, sel = sel, selprops = selprops, selthr = selthr))
}

#a <- Simulation_study(seed = 5, height = 5, n = 1000, pk = 500, ev_xy = 0.1)

#########################
no_cores <- nchunks
cl <- makeCluster(no_cores)

clusterEvalQ(cl, library(rpart))
clusterEvalQ(cl, library(sharp))
clusterEvalQ(cl, library(fake))
clusterEvalQ(cl, library(rlist))
clusterEvalQ(cl, library(randomForest))
clusterEvalQ(cl, library(xgboost))
clusterEvalQ(cl, library(Boruta))
clusterEvalQ(cl, library(rfvimptest))
clusterEvalQ(cl, library(data.tree))
clusterExport(cl, c("numrep", "height", "n", "pk", "ev_xy", "HugeAdjacency",
                    "reorderCV", "CART1", "CART2", "CART3", "getMetrics", "Simulation_study",
                    "SimulateTree", "Bin2Dec", "GenerateTree", "PruneTree", "getImpXgboost_new"))

out <- parLapply(1:numrep,
                 function(i) {Simulation_study(seed = i, height = height, 
                                               n = n, pk = pk, ev_xy = ev_xy)},
                 cl = cl)

stopCluster(cl)

########### cleanup output
Metrics <- sapply(out, "[", 1)
Metrics <- list.rbind(Metrics)
Metrics$model <- factor(Metrics$model, levels = unique(Metrics$model))

Sels <- sapply(out, "[", 2)
Sels <- list.rbind(Sels)

Selprops <- sapply(out, "[", 3)
Selprops <- list.rbind(Selprops)

Selthrs <- sapply(out, "[", 4)
Selthrs <- list.rbind(Selthrs)

theta <- rep(0, pk)
names(theta) <- paste0("var", 1:pk)
numvar <- sum(2^c(0:(height-2)))
theta[1:numvar] <- 1

########### save output
saveRDS(Metrics, file = paste(folder,"Metrics.rds", sep = "/"))
saveRDS(Sels, file = paste(folder,"Sels.rds", sep = "/"))
saveRDS(Selprops, file = paste(folder,"Selprops.rds", sep = "/"))
saveRDS(Selthrs, file = paste(folder,"Selthrs.rds", sep = "/"))
saveRDS(theta, file = paste(folder,"Theta.rds", sep = "/"))

###########
time2 <- Sys.time()
print(difftime(time2, time1), units = "mins")
