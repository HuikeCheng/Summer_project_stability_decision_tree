######## pass in arguments
args=commandArgs(trailingOnly=TRUE)
numrep <- as.numeric(args[1])
n <- as.numeric(args[2])
pk <- as.numeric(args[3])
nu_xy <- as.numeric(args[4])
ev_xy <- as.numeric(args[5])
nchunks <- as.numeric(args[6])

######## create folder
folder <- paste("Linear", n, pk, nu_xy, ev_xy, sep = "_")
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
Simulation_study <- function(seed, n, pk, nu_xy, ev_xy) {
  # set seed
  set.seed(seed)
  # generate data
  simul <- SimulateRegression(n = n, pk = pk, nu_xy = nu_xy, ev_xy = ev_xy,
                              beta_abs = 1, beta_sign = 1)
  # vars
  myvars_T <- rownames(simul$beta)[which(simul$beta != 0)]
  myvars_F <- rownames(simul$beta)[which(simul$beta == 0)]
  
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
  
  ####### Stability-lasso
  lasso <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = PenalisedRegression,
    family = "gaussian")
  
  myvars <- names(SelectedVariables(lasso))[which(SelectedVariables(lasso)!=0)]
  metrics$LASSO <- getMetrics(getCM(vars = myvars, Tvars = myvars_T, Fvars = myvars_F))
  
  ####### output
  return(metrics)
}

# a <- Simulation_study(seed = 1, n = 100, pk = 50, nu_xy = 0.2, ev_xy = 0.5)
# b <- Simulation_study(seed = 2, n = 100, pk = 50, nu_xy = 0.2, ev_xy = 0.5)

########### Parallelise
no_cores <- nchunks
cl <- makeCluster(no_cores)

clusterEvalQ(cl, library(rpart))
clusterEvalQ(cl, library(sharp))
clusterEvalQ(cl, library(fake))
clusterExport(cl, c("numrep", "n", "pk", "nu_xy", "ev_xy",
                    "CART1", "CART2", "getBeta", "getCM", "getMetrics", "Simulation_study"))

out <- pblapply(1:numrep,
                function(i) {Simulation_study(seed = i, n = n, pk = pk, nu_xy = nu_xy, ev_xy = ev_xy)},
                cl = cl)

stopCluster(cl)

########### cleanup output
Metrics <- vector("list", length = numrep)
for (i in seq_along(out)) {
  Metrics[[i]] <- as.data.frame(list.rbind(out[[i]]))
  Metrics[[i]]$model <- rownames(Metrics[[i]])
}
Metrics <- list.rbind(Metrics)
Metrics$model <- factor(Metrics$model, levels = unique(Metrics$model))

########### save output
saveRDS(Metrics, file = paste(folder,"Metrics.rds", sep = "/"))

########### plots - Metrics
title <- paste("Linear: ", "n = ", n, ", pk = ", pk, ", nu_xy = ", nu_xy, ", ev_xy = ", ev_xy, sep = "")

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