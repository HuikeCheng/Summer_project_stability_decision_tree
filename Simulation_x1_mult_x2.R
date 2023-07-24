######## pass in arguments
args=commandArgs(trailingOnly=TRUE)
# args <- c("100", "1000", "500", "0.2", "2","50")
numrep <- as.numeric(args[1])
n <- as.numeric(args[2])
pk <- as.numeric(args[3])
ev_xy <- as.numeric(args[4])
numtrue <- as.numeric(args[5])
nchunks <- as.numeric(args[6])

# job script

######## create folder
folder <- paste(paste0("Outputs/", "X1MultX2"), n, pk, ev_xy, numtrue, sep = "_")
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
library(pbapply)
library(ggplot2)
source("Functions.R")
source("stab1.R")
source("stab2.R")
source("SimulateX.R")

#####################
Simulation_study <- function(seed, n, pk, association, ev_xy, numtrue = 2, multivariate_normal = TRUE) {
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
  if (numtrue == 2) {
    mu <- xdata[,1]*xdata[,2]
  } else if (numtrue > 2) {
    mu <- xdata[,1]
    for (i in 2:numtrue) {
      mu <- mu*xdata[,i]
    }
  } else {
    stop("numtrue must be at least 2")
  }
  
  ypred <- transformX(mu, association)
  sigma <- sqrt((1 - ev_xy) / ev_xy * stats::var(ypred))
  ydata <- stats::rnorm(n = n, mean = ypred, sd = sigma)
  ydata <- as.matrix(ydata)
  
  # theta
  theta <- rep(0, pk)
  names(theta) <- colnames(xdata)
  theta[1:numtrue] <- 1
  
  # set up output
  metrics <- vector("list", length = 5)
  names(metrics) <- c("CART", "CART_ms0", "sCART1", "sCART2", "sLASSO")
  
  ######## run CART
  mydata <- data.frame(ydata = ydata[, 1], xdata)
  cart <- rpart(ydata ~ ., mydata, method = "anova")
  
  myvars <- unique(rownames(cart$splits))
  selvars <- rep(0, pk)
  names(selvars) <- colnames(xdata)
  selvars[which(names(selvars) %in% myvars)] <- 1
  metrics$CART <- getMetrics(selvars = selvars, theta = theta)
  
  ######## run CART-maxsurrogate = 0
  cartms0 <- rpart(ydata ~ ., mydata, method = "anova", maxsurrogate = 0)
  
  myvars <- unique(rownames(cartms0$splits))
  selvars <- rep(0, pk)
  names(selvars) <- colnames(xdata)
  selvars[which(names(selvars) %in% myvars)] <- 1
  metrics$CART_ms0 <- getMetrics(selvars = selvars, theta = theta)
  
  ######### run stab1
  scart1 <- sharp::VariableSelection(
    xdata = xdata, ydata = ydata,
    implementation = CART1,
    family = "gaussian",
    maxsurrogate = 0)
  
  selvars <- SelectedVariables(scart1)
  metrics$sCART1 <- getMetrics(selvars = selvars, theta = theta)
  
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
  
  ####### run Stability-lasso
  slasso <- sharp::VariableSelection(
    xdata = xdata, ydata = ydata,
    implementation = PenalisedRegression,
    family = "gaussian")
  
  selvars <- SelectedVariables(slasso)
  metrics$sLASSO <- getMetrics(selvars = selvars, theta = theta)
  
  ####### Cleanup output
  metrics <- as.data.frame(list.rbind(metrics))
  metrics$model <- rownames(metrics)
  
  ####### output
  return(metrics)
}

# a <- Simulation_study(seed = 1, n = 1000, pk = 500, association = "Linear", ev_xy = 0.2, numtrue = 3)

#########################
Simulation <- function(seed, n, pk, myassoc, ev_xy, numtrue) {
  # set up output
  out <- vector("list", length=length(myassoc))
  # run simulation for each nl
  for (i in 1:length(myassoc)) {
    out[[i]] <- Simulation_study(seed = seed, n = n, pk = pk, association = myassoc[i], 
                                 ev_xy = ev_xy, numtrue = numtrue)
  }
  # clean up
  names(out) <- myassoc
  out <- as.data.frame(list.rbind(out))
  out$relationship <- rep(myassoc, each=5)
  # output
  return(out)
}

#########################
myassoc <- c("Heteroscedastic", "Absolute", "Quadratic", "Cosine", "Sine", "Indicator",
             "SquareWave", "Sawtooth", "Exponential", "Sigmoidal", "Cubic", "Linear")

no_cores <- nchunks
cl <- makeCluster(no_cores)

clusterEvalQ(cl, library(rpart))
clusterEvalQ(cl, library(sharp))
clusterEvalQ(cl, library(fake))
clusterEvalQ(cl, library(rlist))
clusterExport(cl, c("n", "pk", "myassoc", "ev_xy", "numtrue", "numrep",
                    "CART1", "CART2", "getMetrics", "Simulation_study",
                    "HugeAdjacency","transformX", "Simulation"))

out <- pblapply(1:numrep,
                function(i) {Simulation(seed = i, n = n, pk = pk, myassoc = myassoc, 
                                        ev_xy = ev_xy, numtrue = numtrue)},
                cl = cl)

stopCluster(cl)

########### cleanup output
Metrics <- list.rbind(out)
Metrics$model <- factor(Metrics$model, levels = unique(Metrics$model))
Metrics$relationship <- factor(Metrics$relationship, levels = myassoc)

########### save output
saveRDS(Metrics, file = paste(folder,"Metrics.rds", sep = "/"))

########### plots - Metrics
title <- paste("nl(X1*X2): ", "n=", n, ", pk=", pk, ", ev_xy=", ev_xy, ", numtrue=", numtrue, sep = "")

png(file = paste(folder,"Recall.png", sep = "/"))
ggplot(Metrics, aes(x=Recall, fill=forcats::fct_rev(model), y=Relationship)) +
  geom_boxplot() +
  coord_flip() +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom") +
  scale_y_discrete(guide = guide_axis(angle = 45)) +
  ggtitle(title)
dev.off()

png(file = paste(folder,"Precision.png", sep = "/"))
ggplot(Metrics, aes(x=Precision, fill=forcats::fct_rev(model), y=Relationship)) +
  geom_boxplot() +
  coord_flip() +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom") +
  scale_y_discrete(guide = guide_axis(angle = 45)) +
  ggtitle(title)
dev.off()

png(file = paste(folder,"F1.png", sep = "/"))
ggplot(Metrics, aes(x=F1, fill=forcats::fct_rev(model), y=Relationship)) +
  geom_boxplot() +
  coord_flip() +
  labs(x = "F1-score") +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom") +
  scale_y_discrete(guide = guide_axis(angle = 45)) +
  ggtitle(title)
dev.off()

png(file = paste(folder,"Specificity.png", sep = "/"))
ggplot(Metrics, aes(x=Specificity, fill=forcats::fct_rev(model), y=Relationship)) +
  geom_boxplot() +
  coord_flip() +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom") +
  scale_y_discrete(guide = guide_axis(angle = 45)) +
  ggtitle(title)
dev.off()

png(file = paste(folder,"HD.png", sep = "/"))
ggplot(Metrics, aes(x=HD, fill=forcats::fct_rev(model), y=Relationship)) +
  geom_boxplot() +
  coord_flip() +
  labs(x = "Hamming distance") +
theme(legend.title = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "bottom") +
  scale_y_discrete(guide = guide_axis(angle = 45)) +
  ggtitle(title)
dev.off()