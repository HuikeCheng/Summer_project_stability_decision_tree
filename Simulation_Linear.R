######### load packages
library(rpart)
library(sharp)
library(fake)
library(rlist)
library(ggplot2)
source("stab1.R")
source("stab2.R")
source("Functions.R")

####################### Set up simulation study ################################
Simulation_study <- function(seed, n, pk, nu_xy) {
  # set seed
  set.seed(seed)
  # generate data
  simul <- SimulateRegression(n = n, pk = pk, nu_xy = nu_xy)
  # vars
  myvars_T <- rownames(simul$beta)[which(simul$beta != 0)]
  myvars_F <- rownames(simul$beta)[which(simul$beta == 0)]
  
  # set up output
  metrics <- vector("list", length = 8)
  names(metrics) <- c("CART", "stab1", "stab2", "stab1MS1", "stab2MS1", "stab1MS2", "stab2MS2", "LASSO")
  meanBeta <- vector("list", length = 8)
  names(meanBeta) <- c("CART", "stab1", "stab2", "stab1MS1", "stab2MS1", "stab1MS2", "stab2MS2", "LASSO")
  
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
  
  ######### run stab1 with ms = 1
  stab1MS1 <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = CART1,
    family = "gaussian",
    maxsurrogate = 1)

  myvars <- names(SelectedVariables(stab1MS1))[which(SelectedVariables(stab1MS1)!=0)]
  metrics$stab1MS1 <- getMetrics(getCM(vars = myvars, Tvars = myvars_T, Fvars = myvars_F))
  meanBeta$stab1MS1 <- getBeta(vars = myvars, Tvars = myvars_T, beta = simul$beta)
  
  ####### run stab2 with ms = 1
  stab2MS1 <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = CART2,
    family = "gaussian",
    Lambda = Lambda,
    maxsurrogate = 1)

  myvars <- names(SelectedVariables(stab2MS1))[which(SelectedVariables(stab2MS1)!=0)]
  metrics$stab2MS1 <- getMetrics(getCM(vars = myvars, Tvars = myvars_T, Fvars = myvars_F))
  meanBeta$stab2MS1 <- getBeta(vars = myvars, Tvars = myvars_T, beta = simul$beta)
  
  ######### run stab1 with ms = 2
  stab1MS2 <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = CART1,
    family = "gaussian",
    maxsurrogate = 2)

  myvars <- names(SelectedVariables(stab1MS2))[which(SelectedVariables(stab1MS2)!=0)]
  metrics$stab1MS2 <- getMetrics(getCM(vars = myvars, Tvars = myvars_T, Fvars = myvars_F))
  meanBeta$stab1MS2 <- getBeta(vars = myvars, Tvars = myvars_T, beta = simul$beta)
  
  ####### run stab2 with ms = 2
  stab2MS2 <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = CART2,
    family = "gaussian",
    Lambda = Lambda,
    maxsurrogate = 2)

  myvars <- names(SelectedVariables(stab2MS2))[which(SelectedVariables(stab2MS2)!=0)]
  metrics$stab2MS2 <- getMetrics(getCM(vars = myvars, Tvars = myvars_T, Fvars = myvars_F))
  meanBeta$stab2MS2 <- getBeta(vars = myvars, Tvars = myvars_T, beta = simul$beta)
  
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

# a <- Simulation_study(seed = 1, n = 100, pk = 50, nu_xy = 0.2)
# b <- Simulation_study(seed = 2, n = 100, pk = 50, nu_xy = 0.2)

############################## Set up repetitions ##############################
Run_simu <- function(numrep, n, pk, nu_xy) {
  # set up output
  Metrics <- vector("list", length = 8)
  names(Metrics) <- c("CART", "stab1", "stab2", "stab1MS1", "stab2MS1", "stab1MS2", "stab2MS2", "LASSO")
  for (i in seq_along(Metrics)) {
    Metrics[[i]] <- vector("list", length = numrep)
  }
  Betas <- Metrics
  
  # run simulation
  for (i in 1:numrep){
    print(paste0("Simulation", i))
    
    tmp <- Simulation_study(seed = i, n = n, pk = pk, nu_xy = nu_xy)
    
    Metrics$CART[[i]] <- tmp[[1]]$CART
    Metrics$stab1[[i]] <- tmp[[1]]$stab1
    Metrics$stab2[[i]] <- tmp[[1]]$stab2
    Metrics$stab1MS1[[i]] <- tmp[[1]]$stab1MS1
    Metrics$stab2MS1[[i]] <- tmp[[1]]$stab2MS1
    Metrics$stab1MS2[[i]] <- tmp[[1]]$stab1MS2
    Metrics$stab2MS2[[i]] <- tmp[[1]]$stab2MS2
    Metrics$LASSO[[i]] <- tmp[[1]]$LASSO
    
    Betas$CART[[i]] <- tmp[[2]]$CART
    Betas$stab1[[i]] <- tmp[[2]]$stab1
    Betas$stab2[[i]] <- tmp[[2]]$stab2
    Betas$stab1MS1[[i]] <- tmp[[2]]$stab1MS1
    Betas$stab2MS1[[i]] <- tmp[[2]]$stab2MS1
    Betas$stab1MS2[[i]] <- tmp[[2]]$stab1MS2
    Betas$stab2MS2[[i]] <- tmp[[2]]$stab2MS2
    Betas$LASSO[[i]] <- tmp[[2]]$LASSO
  }
  
  # output
  return(list(Metrics = Metrics,
              Betas = Betas))
}

########### Run simulation
out <- Run_simu(numrep = 3, n = 100, pk = 50, nu_xy = 0.2)

########### cleanup output
for (i in seq_along(out$Metrics)) {
  out$Metrics[[i]] <- list.rbind(out$Metrics[[i]])
}

for (i in seq_along(out$Betas)) {
  out$Betas[[i]] <- list.rbind(out$Betas[[i]])
}

########### summary measures
metrics <- vector("list", length = 8)
names(metrics) <- names(out$Metrics)
for (i in seq_along(out$Metrics)) {
  metrics[[i]] <- apply(out$Metrics[[i]], 2, median)
}
metrics <- list.rbind(metrics)
metrics

betas <- vector("list", length = 8)
names(betas) <- names(out$Betas)
for (i in seq_along(out$Betas)) {
  betas[[i]] <- apply(out$Betas[[i]], 2, mean)
}
betas <- list.rbind(betas)
betas

########### plots - Metrics
numrep <- 3
plot_metrics <- as.data.frame(list.rbind(out$Metrics))
plot_metrics$model <- as.factor(rep(names(out$Metrics), each = numrep))
plot_metrics$model <- factor(plot_metrics$model, levels = unique(plot_metrics$model))

ggplot(plot_metrics, aes(x=Recall, y=forcats::fct_rev(model))) +
  geom_boxplot() +
  theme(legend.title = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle("Linear: numrep = 10, n = 100, pk = 50, nu_xy = 0.2")

ggplot(plot_metrics, aes(x=Precision, y=forcats::fct_rev(model))) +
  geom_boxplot() +
  theme(legend.title = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle("Linear: numrep = 10, n = 100, pk = 50, nu_xy = 0.2")

ggplot(plot_metrics, aes(x=F1, y=forcats::fct_rev(model))) +
  geom_boxplot() +
  theme(legend.title = element_blank(),
        axis.title.y = element_blank()) +
  labs(x = "F1-score") +
  ggtitle("Linear: numrep = 10, n = 100, pk = 50, nu_xy = 0.2")

ggplot(plot_metrics, aes(x=Specificity, y=forcats::fct_rev(model))) +
  geom_boxplot() +
  theme(legend.title = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle("Linear: numrep = 10, n = 100, pk = 50, nu_xy = 0.2")


########### plots - Betas
plot_betas <- as.data.frame(list.rbind(out$Betas))
plot_betas$model <- as.factor(rep(names(out$Betas), each = numrep))
plot_betas$model <- factor(plot_betas$model, levels = unique(plot_betas$model))

ggplot(plot_betas, aes(y = meanTP, x= forcats::fct_rev(model))) +
  geom_point() +
  geom_errorbar(aes(ymin=min(minTP), ymax=max(maxTP)), width=.2) +
  ylim(0,1) +
  coord_flip() +
  labs(y = "Beta value of TPs", x = "") +
  ggtitle("Linear: numrep = 10, n = 100, pk = 50, nu_xy = 0.2")

ggplot(plot_betas, aes(y = meanFN, x= forcats::fct_rev(model))) +
  geom_point() +
  geom_errorbar(aes(ymin=min(minFN), ymax=max(maxFN)), width=.2) +
  ylim(0,1) +
  coord_flip() +
  labs(y = "Beta value of FNs", x = "") +
  ggtitle("Linear: numrep = 10, n = 100, pk = 50, nu_xy = 0.2")
