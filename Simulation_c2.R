######## pass in arguments
args=commandArgs(trailingOnly=TRUE)
#args <- c("100", "1000", "10000", "0.01", "5")
numrep <- as.numeric(args[1])
n <- as.numeric(args[2])
pk <- as.numeric(args[3])
ev_xy <- as.numeric(args[4])
nchunks <- as.numeric(args[5])

######## create folder
folder <- paste(paste0("Outputs/", "Complex2"), n, pk, ev_xy, sep = "_")
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
library(pbapply)
library(ggplot2)
source("stab1.R")
source("stab2.R")
source("Functions.R")

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
  
  # theta
  theta <- rep(0, pk)
  names(theta) <- colnames(xdata)
  theta[1:3] <- 1
  
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

# a <- Simulation_study(seed = 1, n = 100, pk = 500, ev_xy = 0.01)

########### Parallelise
no_cores <- nchunks
cl <- makeCluster(no_cores)

clusterEvalQ(cl, library(rpart))
clusterEvalQ(cl, library(sharp))
clusterEvalQ(cl, library(fake))
clusterEvalQ(cl, library(rlist))
clusterExport(cl, c("numrep", "n", "pk", "ev_xy", "HugeAdjacency",
                    "CART1", "CART2", "getMetrics", "Simulation_study"))

out <- pblapply(1:numrep,
                function(i) {Simulation_study(seed = i, n = n, pk = pk, ev_xy = ev_xy)},
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
theta[1:3] <- 1

########### save output
saveRDS(Metrics, file = paste(folder,"Metrics.rds", sep = "/"))
saveRDS(Sels, file = paste(folder,"Sels.rds", sep = "/"))
saveRDS(theta, file = paste(folder,"Theta.rds", sep = "/"))

########### plots - Metrics
title <- paste("Complex2: ", "n = ", n, ", pk = ", pk, ", ev_xy = ", ev_xy, sep = "")

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