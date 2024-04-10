######## pass in arguments
args=commandArgs(trailingOnly=TRUE)
#args <- c("10", "Linear", "100", "100", "0.2", "0.2", "0.2","5")
numrep <- as.numeric(args[1])
association <- as.character(args[2])
n <- as.numeric(args[3])
pk <- as.numeric(args[4])
ev_xy <- as.numeric(args[5])
nu_xy <- as.numeric(args[6])
nu_nl <- as.numeric(args[7])
nchunks <- as.numeric(args[8])

######## create folder
folder <- paste(paste0("Outputs/", "Mult"), association, n, pk, ev_xy, nu_xy, nu_nl, sep = "_")
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
library(randomForest)
library(Boruta)
library(parallel)
library(ggplot2)
library(xgboost)
library(rfvimptest)
source("Functions/Functions.R")
source("Functions/cart1.R")
source("Functions/cart2.R")
source("Functions/SimulateX.R")

#####################
time1 <- Sys.time()

#####################
Simulation_study <- function(seed, n, pk, association, ev_xy, nu_xy, nu_nl, multivariate_normal = TRUE) {
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
  numvar <- nu_xy*pk
  nlvar <- nu_nl*numvar
  lvar <- numvar - nlvar
  
  if (nu_nl == 1) {
    mu <- xdata[,51]
    for (i in 2:numvar) {
      mu <- mu*xdata[,(50+i)]
    }
    ypred <- transformX(mu, association)
  } else {
    beta_l <- rep(1, lvar)
    mu_l <- as.matrix(xdata[,51:(50+lvar)])%*%beta_l
    
    mu_nl <- xdata[,(51+lvar)]
    for (i in 2:nlvar) {
      mu_nl <- mu_nl*xdata[,(50+lvar+i)]
    }
    ypred <- mu_l + transformX(mu_nl, association)
  }
  
  sigma <- sqrt((1 - ev_xy) / ev_xy * stats::var(ypred))
  ydata <- stats::rnorm(n = n, mean = ypred, sd = sigma)
  ydata <- as.matrix(ydata)
  
  # theta
  theta <- rep(0, pk)
  names(theta) <- colnames(xdata)
  theta[51:(50+numvar)] <- 1
  
  # set up output
  metrics <- vector("list", length = 10)
  sel <- vector("list", length = 10)
  names(metrics) <- c("CART", "sCART1", "sCART2", "sLASSO", 
                      "RF_B_CT", "RF_B_C", "XGB_B_CT", "XGB_B_C", 
                      "RF_VT_SRPT", "RF_VT_SAPT")
  names(sel) <- c("CART", "sCART1", "sCART2", "sLASSO", 
                  "RF_B_CT", "RF_B_C", "XGB_B_CT", "XGB_B_C", 
                  "RF_VT_SRPT", "RF_VT_SAPT")
  
  selprops <- vector("list", length = 3)
  selthr <- vector("list", length = 3)
  names(selprops) <- c("sCART1", "sCART2", "sLASSO")
  names(selthr) <- c("sCART1", "sCART2", "sLASSO")
  
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
    maxsurrogate = 0)
  
  selvars <- SelectedVariables(scart1)
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
  metrics$sCART2 <- getMetrics(selvars = selvars, theta = theta)
  sel$sCART2 <- selvars
  
  selprops$sCART2 <- scart2$selprop[ArgmaxId(scart2)[1],]
  selthr$sCART2 <- Argmax(scart2)[2]
  
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
  
  ####### run RF_VT
  # rf_vt_sprt <- rfvimptest_new(data = mydata, yname = "ydata", test = "general", 
  #                              type = "SPRT")
  # myvars <- names(rf_vt_sprt$testres)[which(rf_vt_sprt$testres == "accept H1")]
  # selvars <- rep(0, pk)
  # names(selvars) <- colnames(xdata)
  # selvars[which(names(selvars) %in% myvars)] <- 1
  # metrics$RF_VT_SPRT <- getMetrics(selvars = selvars, theta = theta)
  # sel$RF_VT_SPRT <- selvars
  # 
  # rf_vt_sapt <- rfvimptest_new(data = mydata, yname = "ydata", test = "general", 
  #                              type = "SAPT")
  # myvars <- names(rf_vt_sapt$testres)[which(rf_vt_sapt$testres == "accept H1")]
  # selvars <- rep(0, pk)
  # names(selvars) <- colnames(xdata)
  # selvars[which(names(selvars) %in% myvars)] <- 1
  # metrics$RF_VT_SAPT <- getMetrics(selvars = selvars, theta = theta)
  # sel$RF_VT_SAPT <- selvars
  
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

#a <- Simulation_study(seed = 1, n = 100, pk = 100, association = "Exponential", ev_xy = 0.2, nu_xy = 0.2, nu_nl = 0.2)

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
clusterExport(cl, c("n", "pk", "association","ev_xy", "nu_xy", "nu_nl", "numrep",
                    "CART1", "CART2", "getMetrics", "Simulation_study",
                    "HugeAdjacency","transformX", "getImpXgboost_new"))

out <- parLapply(1:numrep,
                 function(i) {Simulation_study(seed = i, n = n, pk = pk, association = association,
                                               ev_xy = ev_xy, nu_xy = nu_xy, nu_nl = nu_nl)},
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
theta[51:(50+nu_xy*pk)] <- 1

########### save output
saveRDS(Metrics, file = paste(folder,"Metrics.rds", sep = "/"))
saveRDS(Sels, file = paste(folder,"Sels.rds", sep = "/"))
saveRDS(Selprops, file = paste(folder,"Selprops.rds", sep = "/"))
saveRDS(Selthrs, file = paste(folder,"Selthrs.rds", sep = "/"))
saveRDS(theta, file = paste(folder,"Theta.rds", sep = "/"))

###########
time2 <- Sys.time()
print(difftime(time2, time1), units = "mins")