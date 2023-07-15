######## pass in arguments
args=commandArgs(trailingOnly=TRUE)
#args <- c("100", "1000", "500", "0.1", "5")
numrep <- as.numeric(args[1])
n <- as.numeric(args[2])
pk <- as.numeric(args[3])
sd <- as.numeric(args[4])
nchunks <- as.numeric(args[5])

######## create folder
folder <- paste(paste0("Outputs/", Complex1), n, pk, sd, sep = "_")
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
library(parallel)
library(pbapply)
library(ggplot2)
source("stab1.R")
source("stab2.R")
source("Functions.R")

####################### Set up simulation study ################################
Simulation_study <- function(seed, n, pk, association, sd) {
  # set seed
  set.seed(seed)
  # simulate data
  X <- matrix(rnorm(n*pk),n,pk)
  Y <- 0.5*X[,1] + 0.45*X[,2] + 0.3*X[,3] + 
    1.5*ifelse(X[,1]<=0 & X[,2]>0 & X[,3] <= 0, 1, 0) +
    0.25*ifelse(X[,1]<=0 & X[,3]>0, 1, 0) + 
    0.14*ifelse(X[,1]>0 & X[,2] > 0, 1, 0) + 
    rnorm(n,0,sd)
  colnames(X) <- paste0("var", 1:pk)
  simul <- list(ydata=as.matrix(Y), xdata=X)
  
  # vars
  myvars_T <- paste0("var", c(1,2,3))
  myvars_F <- paste0("var", c(4:pk))
  
  # set up output
  metrics <- vector("list", length = 4)
  names(metrics) <- c("CART", "stab1", "stab2", "LASSO")
  
  ######## run CART
  mydata <- data.frame(ydata = simul$ydata[, 1], simul$xdata)
  mycart <- rpart(ydata ~ ., mydata, method="anova")
  
  myvars <- unique(rownames(mycart$splits))
  metrics$CART <- getMetrics(getCM(vars = myvars, Tvars = myvars_T, Fvars = myvars_F))
  
  ######### run stab1
  stab1 <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = CART1,
    family = "gaussian",
    maxsurrogate = 0)
  
  myvars <- names(SelectedVariables(stab1))[which(SelectedVariables(stab1)!=0)]
  metrics$stab1 <- getMetrics(getCM(vars = myvars, Tvars = myvars_T, Fvars = myvars_F))
  
  ####### run stab2
  Lambda <- mycart$cptable[,1]
  stab2 <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = CART2,
    family = "gaussian",
    Lambda = Lambda,
    maxsurrogate = 0)
  
  myvars <- names(SelectedVariables(stab2))[which(SelectedVariables(stab2)!=0)]
  metrics$stab2 <- getMetrics(getCM(vars = myvars, Tvars = myvars_T, Fvars = myvars_F))
  
  ####### run Stability-lasso
  lasso <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = PenalisedRegression,
    family = "gaussian")
  
  myvars <- names(SelectedVariables(lasso))[which(SelectedVariables(lasso)!=0)]
  metrics$LASSO <- getMetrics(getCM(vars = myvars, Tvars = myvars_T, Fvars = myvars_F))
  
  ####### Cleanup output
  metrics <- as.data.frame(list.rbind(metrics))
  metrics$model <- rownames(metrics)
  
  ####### output
  return(metrics)
}

# a <- Simulation_study(seed = 1, n = 100, pk = 50, sd = 0.1)

########### Parallelise
no_cores <- nchunks
cl <- makeCluster(no_cores)

clusterEvalQ(cl, library(rpart))
clusterEvalQ(cl, library(sharp))
clusterEvalQ(cl, library(fake))
clusterExport(cl, c("numrep", "n", "pk", "sd",
                    "CART1", "CART2", "getBeta", "getCM", "getMetrics", "Simulation_study"))

out <- pblapply(1:numrep,
                function(i) {Simulation_study(seed = i, association = association, 
                                              n = n, pk = pk, sd = sd)},
                cl = cl)

stopCluster(cl)

########### cleanup output
Metrics <- list.rbind(Metrics)
Metrics$model <- factor(Metrics$model, levels = unique(Metrics$model))

########### save output
saveRDS(Metrics, file = paste(folder,"Metrics.rds", sep = "/"))

########### plots - Metrics
title <- paste("Complex1: ", "n = ", n, ", pk = ", pk, ", sd = ", sd, sep = "")

png(file = paste(folder,"Recall.png", sep = "/"))
ggplot(Metrics, aes(x=Recall, y=forcats::fct_rev(model))) +
  geom_boxplot() +
  theme(legend.title = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle(title)
dev.off()

png(file = paste(folder,"Precision.png", sep = "/"))
ggplot(Metrics, aes(x=Precision, y=forcats::fct_rev(model))) +
  geom_boxplot() +
  theme(legend.title = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle(title)
dev.off()

png(file = paste(folder,"F1.png", sep = "/"))
ggplot(Metrics, aes(x=F1, y=forcats::fct_rev(model))) +
  geom_boxplot() +
  theme(legend.title = element_blank(),
        axis.title.y = element_blank()) +
  labs(x = "F1-score") +
  ggtitle(title)
dev.off()

png(file = paste(folder,"Specificity.png", sep = "/"))
ggplot(Metrics, aes(x=Specificity, y=forcats::fct_rev(model))) +
  geom_boxplot() +
  theme(legend.title = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle(title)
dev.off()