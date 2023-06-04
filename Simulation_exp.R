######### load packages
library(rpart)
library(sharp)
library(fake)
library(rlist)
library(tidyverse)
library(ggplot2)
source("stab1.R")
source("stab2.R")
source("Functions.R")

####################### Set up simulation study ################################
Simulation_study <- function(seed, n, pk, nu_xy, proportion) {
  # set seed
  set.seed(seed)
  # generate data
  simul <- SimulateRegression(n = n, pk = pk, nu_xy = nu_xy)
  # vars
  myvars_T <- rownames(simul$beta)[which(simul$beta != 0)]
  myvars_F <- rownames(simul$beta)[which(simul$beta == 0)]
  # get the true vars
  id <- which(colnames(simul$xdata) %in% myvars_T)
  id <- sample(id, round(length(id)*proportion, 0))
  # change these vars to abs(X)
  for (i in 1:length(id)){
    simul$xdata[,id[i]] <- abs(simul$xdata[,id[i]])
  }
  simul$ydata <- simul$xdata%*%simul$beta
  # change these vars to log(X)
  for (i in 1:length(id)){
    simul$xdata[,id[i]] <- log(simul$xdata[,id[i]])
  }
  # set up output
  metrics <- vector("list", length = 4)
  names(metrics) <- c("CART", "stab1", "stab2", "LASSO")
  meanBeta <- vector("list", length = 4)
  names(meanBeta) <- c("CART", "stab1", "stab2", "LASSO")
  
  ######## run CART
  mydata <- data.frame(ydata = simul$ydata[, 1], simul$xdata)
  mycart <- rpart(ydata ~ ., mydata, method="anova")
  
  myvars <- names(mycart$variable.importance)
  metrics$CART <- getMetrics(getCM(vars = myvars, Tvars = myvars_T, Fvars = myvars_F))
  meanBeta$CART <- getBeta(vars = myvars, Tvars = myvars_T, beta = simul$beta)
  
  ######### run stab1
  stab1 <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = CART1,
    family = "gaussian",
    maxsurrogate = 0)
  
  myvars <- names(SelectedVariables(stab1))[which(SelectedVariables(stab1)!=0)]
  metrics$stab1 <- getMetrics(getCM(vars = myvars, Tvars = myvars_T, Fvars = myvars_F))
  meanBeta$stab1 <- getBeta(vars = myvars, Tvars = myvars_T, beta = simul$beta)
  
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
  meanBeta$stab2 <- getBeta(vars = myvars, Tvars = myvars_T, beta = simul$beta)
  
  ####### Stability-lasso
  lasso <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = PenalisedRegression,
    family = "gaussian")
  
  myvars <- names(SelectedVariables(lasso))[which(SelectedVariables(lasso)!=0)]
  metrics$LASSO <- getMetrics(getCM(vars = myvars, Tvars = myvars_T, Fvars = myvars_F))
  meanBeta$LASSO <- getBeta(vars = myvars, Tvars = myvars_T, beta = simul$beta)
  
  ####### output
  return(list(metrics = metrics,
              meanBeta = meanBeta))
}

# a <- Simulation_study(seed = 1, n = 100, pk = 50, nu_xy = 0.2, proportion = 0.5)

########### Parallelise
library(foreach)
library(doParallel)

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

clusterEvalQ(cl, library(fake))
clusterEvalQ(cl, library(rpart))
clusterEvalQ(cl, library(sharp))
clusterExport(cl, c("CART1", "CART2", "getBeta", "getCM", "getMetrics", "Simulation_study"))

out <- foreach(seed = 1:1000, 
               .combine = list,
               .multicombine = TRUE,
               .verbose = TRUE)  %dopar%  
  Simulation_study(seed = seed, n = 1000, pk = 500, nu_xy = 0.1, proportion = 0.5)

stopImplicitCluster()

########### cleanup output
numrep <- 1000
Metrics <- vector("list", length = numrep)
for (i in seq_along(out)) {
  Metrics[[i]] <- as.data.frame(list.rbind(out[[i]]$metrics))
  Metrics[[i]]$model <- rownames(Metrics[[i]])
}
Metrics <- list.rbind(Metrics)

Betas <- vector("list", length = numrep)
for (i in seq_along(out)) {
  Betas[[i]] <- as.data.frame(list.rbind(out[[i]]$meanBeta))
  Betas[[i]]$model <- rownames(Betas[[i]])
}
Betas <- list.rbind(Betas)

########### summary measures
Metrics %>% group_by(model) %>% summarise_all(median)

Betas %>% group_by(model) %>% summarise_all(mean)

########### plots - Metrics
Metrics$model <- factor(Metrics$model, levels = unique(Metrics$model))
title <- "Exp: numrep = 1000, n = 1000, pk = 500, nu_xy = 0.1, proportion = 0.5"

ggplot(Metrics, aes(x=Recall, y=forcats::fct_rev(model))) +
  geom_boxplot() +
  theme(legend.title = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle(title)

ggplot(Metrics, aes(x=Precision, y=forcats::fct_rev(model))) +
  geom_boxplot() +
  theme(legend.title = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle(title)

ggplot(Metrics, aes(x=F1, y=forcats::fct_rev(model))) +
  geom_boxplot() +
  theme(legend.title = element_blank(),
        axis.title.y = element_blank()) +
  labs(x = "F1-score") +
  ggtitle(title)

ggplot(Metrics, aes(x=Specificity, y=forcats::fct_rev(model))) +
  geom_boxplot() +
  theme(legend.title = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle(title)

########### plots - Betas
Betas$model <- factor(Betas$model, levels = unique(Betas$model))

ggplot(Betas, aes(y = meanTP, x= forcats::fct_rev(model))) +
  geom_point() +
  geom_errorbar(aes(ymin=min(minTP), ymax=max(maxTP)), width=.2) +
  ylim(0,1) +
  coord_flip() +
  labs(y = "Beta value of TPs", x = "") +
  ggtitle(title)

ggplot(Betas, aes(y = meanFN, x= forcats::fct_rev(model))) +
  geom_point() +
  geom_errorbar(aes(ymin=min(minFN), ymax=max(maxFN)), width=.2) +
  ylim(0,1) +
  coord_flip() +
  labs(y = "Beta value of FNs", x = "") +
  ggtitle(title)
