######### load packages
library(rpart)
library(sharp)
library(fake)
library(rlist)
library(ggplot2)
source("stab2.R")

######### function
Simulation_study <- function(seed, pk) {
  # set seed
  set.seed(seed)
  # simulate dataset
  simul <- SimulateRegression(pk = pk)
  # vars
  myvars_T <- rownames(simul$beta)[which(simul$beta != 0)]
  myvars_F <- rownames(simul$beta)[which(simul$beta == 0)]
  
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
  
  sensitivity_stab <- TP_stab/length(myvars_T)
  precision_stab <- TP_stab/length(myvars_stab)
  f1_stab <- 2*sensitivity_stab*precision_stab/(sensitivity_stab+precision_stab)
  
  sensitivity_stab2 <- TP_stab2/length(myvars_T)
  precision_stab2 <- TP_stab2/length(myvars_stab2)
  f1_stab2 <- 2*sensitivity_stab2*precision_stab2/(sensitivity_stab2+precision_stab2)
  
  ##### output
  return(list(cart = list(Metrics = metrics_cart,
                        Scores = myscore_cart),
              stab = list(Metrics = metrics_stab,
                          Scores = myscore_stab),
              stab2 = list(Metrics = metrics_stab2,
                           Scores = myscore_stab2)))
}

a <- Simulation_study(1,50)

######### Repeat the simulation 100 times
# set up output
Metrics_stab <- vector("list", length=100)
Metrics_cart <- vector("list", length=100)
Scores_stab <- vector("list", length=100)
Scores_cart <- vector("list", length=100)

# run simulation
for (i in 1:100){
  print(paste0("Simulation", i))
  tmp <- Simulation_study(i,50)
  Metrics_stab[[i]] <- tmp$stab$Metrics
  Metrics_cart[[i]] <- tmp$cart$Metrics
  Scores_stab[[i]] <- tmp$stab$Scores
  Scores_cart[[i]] <- tmp$cart$Scores
}

########### cleanup output
Stab_metrics <- list.rbind(Metrics_stab)
CART_metrics <- list.rbind(Metrics_cart)
summary(Stab_metrics)
summary(CART_metrics)

########### plots
df_metrics <- as.data.frame(rbind(Stab_metrics,CART_metrics))
df_metrics$model <- as.factor(c(rep("Stability", 100), rep("CART", 100)))

ggplot(df_metrics, aes(x=Sensitivity, group=model)) +
  geom_boxplot(aes(fill=model)) +
  theme(legend.title = element_blank(), 
        axis.text.y = element_blank())

ggplot(df_metrics, aes(x=Precision, group=model)) +
  geom_boxplot(aes(fill=model)) +
  theme(legend.title = element_blank(),
        axis.text.y = element_blank())

ggplot(df_metrics, aes(x=F1, group=model)) +
  geom_boxplot(aes(fill=model)) +
  theme(legend.title = element_blank(),
        axis.text.y = element_blank()) +
  labs(x="F1-score")
