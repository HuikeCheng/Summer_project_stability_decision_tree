######### load packages
library(rpart)
library(sharp)
library(fake)
library(rlist)
library(ggplot2)
source("stab2.R")

######### function
Simulation_study <- function(seed, pk, proportion) {
  # set seed
  set.seed(seed)
  # simulate dataset
  simul <- SimulateRegression(pk = pk)
  # vars
  myvars_T <- rownames(simul$beta)[which(simul$beta != 0)]
  myvars_F <- rownames(simul$beta)[which(simul$beta == 0)]
  # get the true vars
  id <- which(colnames(simul$xdata) %in% myvars_T)
  id <- sample(id, round(length(id)*proportion, 0))
  # change these vars to exp(X)
  for (i in 1:length(id)){
    simul$xdata[,id[i]] <- (simul$xdata[,id[i]])^2
  }
  ######## run CART
  mydata <- data.frame(ydata = simul$ydata[, 1], simul$xdata)
  mycart <- rpart(ydata ~ ., mydata, method="anova", )
  myvars_cart <- names(mycart$variable.importance)
  myscore_cart <- mycart$variable.importance
  
  ######### run stab
  stab <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = CART,
    family = "gaussian"
  )
  myvars_stab <- names(SelectedVariables(stab))[which(SelectedVariables(stab)!=0)]
  myscore_stab <- SelectionProportions(stab)
  
  ####### run stab2
  Lambda <- mycart$cptable[,1]
  stab2 <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = CART2,
    family = "gaussian",
    Lambda = Lambda
  )
  myvars_stab2 <- names(SelectedVariables(stab2))[which(SelectedVariables(stab2)!=0)]
  myscore_stab2 <- SelectionProportions(stab2)
  
  ####### Confusion matrix
  TP_cart <- sum(myvars_T %in% myvars_cart)
  FP_cart <- sum(myvars_F %in% myvars_cart)
  TN_cart <- sum(!(myvars_F %in% myvars_cart))
  FN_cart <- sum(!(myvars_T %in% myvars_cart))
  
  TP_stab <- sum(myvars_T %in% myvars_stab)
  FP_stab <- sum(myvars_F %in% myvars_stab)
  TN_stab <- sum(!(myvars_F %in% myvars_stab))
  FN_stab <- sum(!(myvars_T %in% myvars_stab))
  
  TP_stab2 <- sum(myvars_T %in% myvars_stab2)
  FP_stab2 <- sum(myvars_F %in% myvars_stab2)
  TN_stab2 <- sum(!(myvars_F %in% myvars_stab2))
  FN_stab2 <- sum(!(myvars_T %in% myvars_stab2))
  
  ########## metrics
  sensitivity_cart <- TP_cart/length(myvars_T)
  precision_cart <- TP_cart/length(myvars_cart)
  f1_cart <- 2*sensitivity_cart*precision_cart/(sensitivity_cart+precision_cart)
  metrics_cart <- c(sensitivity_cart, precision_cart, f1_cart)
  names(metrics_cart) <- c("Sensitivity", "Precision", "F1")
  
  sensitivity_stab <- TP_stab/length(myvars_T)
  precision_stab <- TP_stab/length(myvars_stab)
  f1_stab <- 2*sensitivity_stab*precision_stab/(sensitivity_stab+precision_stab)
  metrics_stab <- c(sensitivity_stab, precision_stab, f1_stab)
  names(metrics_stab) <- c("Sensitivity", "Precision", "F1")
  
  sensitivity_stab2 <- TP_stab2/length(myvars_T)
  precision_stab2 <- TP_stab2/length(myvars_stab2)
  f1_stab2 <- 2*sensitivity_stab2*precision_stab2/(sensitivity_stab2+precision_stab2)
  metrics_stab2 <- c(sensitivity_stab2, precision_stab2, f1_stab2)
  names(metrics_stab2) <- c("Sensitivity", "Precision", "F1")
  
  ##### results
  Result <- list(cart = list(Metrics = metrics_cart,
                             Scores = myscore_cart),
                 stab = list(Metrics = metrics_stab,
                             Scores = myscore_stab),
                 stab2 = list(Metrics = metrics_stab2,
                              Scores = myscore_stab2))
  
  ########### S vs F1 for stab
  F1 <- matrix(0, nrow=nrow(stab$S_2d), ncol=ncol(stab$S_2d))
  for (i in 1:nrow(stab$S_2d)) {
    for (j in 1:ncol(stab$S_2d)) {
      myvars_stab <- Stable(stab, argmax_id = c(i,j))
      myvars_stab <- names(myvars_stab)[myvars_stab == 1]
      TP_stab <- sum(myvars_T %in% myvars_stab)
      FP_stab <- sum(myvars_F %in% myvars_stab)
      TN_stab <- sum(!(myvars_F %in% myvars_stab))
      FN_stab <- sum(!(myvars_T %in% myvars_stab))
      sensitivity_stab <- TP_stab/length(myvars_T)
      precision_stab <- TP_stab/length(myvars_stab)
      f1_stab <- 2*sensitivity_stab*precision_stab/(sensitivity_stab+precision_stab)
      F1[i,j] <- f1_stab
    }
  }
  
  F1 <- c(F1)
  S <- c(stab$S_2d)
  df <- data.frame(S = S, F1 = F1)
  df <- df[!is.nan(df$S),]
  
  p1 <- ggplot(df, aes(x=F1, y=S)) + geom_point() + ggtitle("stab")
  
  ########### S vs F1 for stab2
  F1 <- matrix(0, nrow=nrow(stab2$S_2d), ncol=ncol(stab2$S_2d))
  for (i in 1:nrow(stab2$S_2d)) {
    for (j in 1:ncol(stab2$S_2d)) {
      myvars_stab2 <- Stable(stab2, argmax_id = c(i,j))
      myvars_stab2 <- names(myvars_stab2)[myvars_stab2 == 1]
      TP_stab2 <- sum(myvars_T %in% myvars_stab2)
      FP_stab2 <- sum(myvars_F %in% myvars_stab2)
      TN_stab2 <- sum(!(myvars_F %in% myvars_stab2))
      FN_stab2 <- sum(!(myvars_T %in% myvars_stab2))
      sensitivity_stab2 <- TP_stab2/length(myvars_T)
      precision_stab2 <- TP_stab2/length(myvars_stab2)
      f1_stab2 <- 2*sensitivity_stab2*precision_stab2/(sensitivity_stab2+precision_stab2)
      F1[i,j] <- f1_stab2
    }
  }
  
  F1 <- c(F1)
  S <- c(stab2$S_2d)
  df <- data.frame(S = S, F1 = F1)
  df <- df[!is.nan(df$S),]
  
  p2 <- ggplot(df, aes(x=F1, y=S)) + geom_point() + ggtitle("stab2")
  
  ####### output
  return(list(Result = Result, stab_SvsF1 = p1, stab2_SvsF1 = p2))
}

a <- Simulation_study(seed = 1, pk = 50, proportion = 1)
a$stab2_SvsF1
a$stab_SvsF1

######### Repeat the simulation 10 times
Run_simu <- function(n, pk, proportion) {
  # set up output
  Metrics_cart <- vector("list", length=n)
  Metrics_stab <- vector("list", length=n)
  Metrics_stab2 <- vector("list", length=n)
  Scores_cart <- vector("list", length=n)
  Scores_stab <- vector("list", length=n)
  Scores_stab2 <- vector("list", length=n)
  SvsF1_stab <- vector("list", length=n)
  SvsF1_stab2 <- vector("list", length=n)
  
  # run simulation
  for (i in 1:n){
    print(paste0("Simulation", i))
    tmp <- Simulation_study(seed = i,pk = pk, proportion = proportion)
    Metrics_cart[[i]] <- tmp$Result$cart$Metrics
    Metrics_stab[[i]] <- tmp$Result$stab$Metrics
    Metrics_stab2[[i]] <- tmp$Result$stab2$Metrics
    Scores_cart[[i]] <- tmp$Result$cart$Scores
    Scores_stab[[i]] <- tmp$Result$stab$Scores
    Scores_stab2[[i]] <- tmp$Result$stab2$Scores
    SvsF1_stab[[i]] <- tmp$stab_SvsF1
    SvsF1_stab2[[i]] <- tmp$stab2_SvsF1
  }
  # output
  return(list(Metrics_cart = Metrics_cart,
              Metrics_stab = Metrics_stab,
              Metrics_stab2 = Metrics_stab2,
              Scores_cart = Scores_cart,
              Scores_stab = Scores_stab,
              Scores_stab2 = Scores_stab2,
              SvsF1_stab = SvsF1_stab,
              SvsF1_stab2 = SvsF1_stab2))
}

########### run
out <- Run_simu(n = 10, pk = 50, proportion = 1)
########### cleanup output
CART_metrics <- list.rbind(out$Metrics_cart)
Stab_metrics <- list.rbind(out$Metrics_stab)
Stab2_metrics <- list.rbind(out$Metrics_stab2)

summary(CART_metrics)
summary(Stab_metrics)
summary(Stab2_metrics)

########### plots
df_metrics <- as.data.frame(rbind(CART_metrics, Stab_metrics, Stab2_metrics))
df_metrics$model <- as.factor(c(rep("CART", 10), rep("Stab", 10), rep("Stab2", 10)))

ggplot(df_metrics, aes(x=Sensitivity, group=model)) +
  geom_boxplot(aes(fill=model)) +
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank()) +
  ggtitle("sr-50-10-0.5")

ggplot(df_metrics, aes(x=Precision, group=model)) +
  geom_boxplot(aes(fill=model)) +
  theme(legend.title = element_blank(),
        axis.text.y = element_blank()) +
  ggtitle("sr-50-10-0.5")

ggplot(df_metrics, aes(x=F1, group=model)) +
  geom_boxplot(aes(fill=model)) +
  theme(legend.title = element_blank(),
        axis.text.y = element_blank()) +
  labs(x="F1-score") +
  ggtitle("sr-50-10-0.5")

######### SvsF1
out$SvsF1_stab[[1]]
out$SvsF1_stab2[[1]]
