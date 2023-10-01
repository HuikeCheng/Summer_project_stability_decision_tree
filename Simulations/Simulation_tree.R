######## pass in arguments
args=commandArgs(trailingOnly=TRUE)
#args <- c("5", "5", "500", "1000", "0.1", "5")
numrep <- as.numeric(args[1])
height <- as.numeric(args[2])
n <- as.numeric(args[3])
pk <- as.numeric(args[4])
ev_xy <- as.numeric(args[5])
nchunks <- as.numeric(args[6])

######## create folder
folder <- paste(paste0("Outputs/Tree"), height, n, pk, ev_xy, sep = "_")
if (file.exists(folder) == FALSE) {
  dir.create(folder)
}
print(folder)

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
source("SimulateTree.R")

####################### Set up simulation study ################################
Simulation_study <- function(seed, height, n, pk, ev_xy) {
  # set seed
  set.seed(seed)
  # generate data
  simul <- SimulateTree(height = height, n = n, pk = pk, ev_xy = ev_xy)
  # theta
  theta <- simul$theta
  
  # set up output
  metrics <- vector("list", length = 6)
  sel <- vector("list", length = 6)
  names(metrics) <- c("CART", "sCART1", "sCART1t","sCART2", "sCART2t","sLASSO")
  names(sel) <- c("CART", "sCART1", "sCART1t","sCART2", "sCART2t","sLASSO")
  
  ######## run CART
  mydata <- data.frame(ydata = simul$ydata[, 1], simul$xdata)
  cart <- rpart(ydata ~ ., mydata, method = "anova", maxsurrogate = 0)
  
  myvars <- unique(rownames(cart$splits))
  selvars <- rep(0, pk)
  names(selvars) <- colnames(simul$xdata)
  selvars[which(names(selvars) %in% myvars)] <- 1
  metrics$CART <- getMetrics(selvars = selvars, theta = theta)
  sel$CART <- selvars
  
  ######### run stab1
  scart1 <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = CART1,
    family = "gaussian",
    maxsurrogate = 0)
  
  selvars <- SelectedVariables(scart1)
  metrics$sCART1 <- getMetrics(selvars = selvars, theta = theta)
  sel$sCART1 <- selvars
  
  ######### run stab1_k=1000
  scart1 <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = CART1,
    family = "gaussian",
    maxsurrogate = 0,
    K=1000)
  
  selvars <- SelectedVariables(scart1)
  metrics$sCART1t <- getMetrics(selvars = selvars, theta = theta)
  sel$sCART1t <- selvars
  
  ####### run stab2
  Lambda <- cart$cptable[,1]
  scart2 <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = CART2,
    family = "gaussian",
    Lambda = Lambda,
    maxsurrogate = 0)
  
  selvars <- SelectedVariables(scart2)
  metrics$sCART2 <- getMetrics(selvars = selvars, theta = theta)
  sel$sCART2 <- selvars
  
  ####### run stab2, k=1000
  scart2 <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = CART2,
    family = "gaussian",
    Lambda = Lambda,
    maxsurrogate = 0,
    K=1000)
  
  selvars <- SelectedVariables(scart2)
  metrics$sCART2t <- getMetrics(selvars = selvars, theta = theta)
  sel$sCART2t <- selvars
  
  ####### run Stability-lasso
  slasso <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
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

a <- Simulation_study(seed = 5, height = 5, n = 1000, pk = 500, ev_xy = 0.1)

########### Parallelise
no_cores <- nchunks
cl <- makeCluster(no_cores)

clusterEvalQ(cl, library(rpart))
clusterEvalQ(cl, library(sharp))
clusterEvalQ(cl, library(fake))
clusterEvalQ(cl, library(rlist))
clusterEvalQ(cl, library(data.tree))
clusterExport(cl, c("numrep", "height", "n", "pk", "ev_xy", "HugeAdjacency",
                    "reorderCV", "CART1", "CART2", "getMetrics", "Simulation_study",
                    "SimulateTree", "Bin2Dec"))

out <- pblapply(1:numrep,
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

########### save output
saveRDS(Metrics, file = paste(folder,"Metrics.rds", sep = "/"))
saveRDS(Sels, file = paste(folder,"Sels.rds", sep = "/"))

########### plots - Metrics
title <- paste("Tree: ", "height = ", height, ", n = ", n, ", pk = ", pk, ", ev_xy = ", ev_xy, sep = "")

png(file = paste(folder,"Recall.png", sep = "/"))
ggplot(Metrics, aes(x=Recall, y=model)) +
  geom_boxplot() +
  coord_flip() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 11)) +
  ggtitle(title)
dev.off()

png(file = paste(folder,"Precision.png", sep = "/"))
ggplot(Metrics, aes(x=Precision, y=model)) +
  geom_boxplot() +
  coord_flip() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 11)) +
  ggtitle(title)
dev.off()

png(file = paste(folder,"F1.png", sep = "/"))
ggplot(Metrics, aes(x=F1, y=model)) +
  geom_boxplot() +
  coord_flip() +
  labs(x = "F1-score") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 11)) +
  ggtitle(title)
dev.off()

png(file = paste(folder,"Specificity.png", sep = "/"))
ggplot(Metrics, aes(x=Specificity, y=model)) +
  geom_boxplot() +
  coord_flip() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 11)) +
  ggtitle(title)
dev.off()

png(file = paste(folder,"HD.png", sep = "/"))
ggplot(Metrics, aes(x=HD, y=model)) +
  geom_boxplot() +
  coord_flip() +
  labs(x = "Hamming distance") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 11)) +
  ggtitle(title)
dev.off()

png(file = paste(folder,"MCC.png", sep = "/"))
ggplot(Metrics, aes(x=MCC, y=model)) +
  geom_boxplot() +
  coord_flip() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 11)) +
  ggtitle(title)
dev.off()