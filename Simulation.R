######### load packages
library(rpart)
library(sharp)
library(fake)
library(rlist)
source("stab1.R")
source("stab2.R")
source("Functions.R")

######### generate data
set.seed(1)
simul <- SimulateRegression(n = 100, pk = 50)

######### true vars
myvars_T <- rownames(simul$beta)[which(simul$beta != 0)]
myvars_F <- rownames(simul$beta)[which(simul$beta == 0)]

######## run CART
mydata <- data.frame(ydata = simul$ydata[, 1], simul$xdata)
mycart <- rpart(ydata ~ ., mydata, method="anova")
myvars_CART <- names(mycart$variable.importance)
# myscore_cart <- mycart$variable.importance

######### run stab1
stab1 <- sharp::VariableSelection(
  xdata = simul$xdata, ydata = simul$ydata,
  implementation = CART1,
  family = "gaussian",
  maxsurrogate = 0)

myvars_stab1 <- names(SelectedVariables(stab1))[which(SelectedVariables(stab1)!=0)]
# myscore_stab1 <- SelectionProportions(stab1)

####### run stab2
Lambda <- mycart$cptable[,1]
stab2 <- sharp::VariableSelection(
  xdata = simul$xdata, ydata = simul$ydata,
  implementation = CART2,
  family = "gaussian",
  Lambda = Lambda,
  maxsurrogate = 0)

myvars_stab2 <- names(SelectedVariables(stab2))[which(SelectedVariables(stab2)!=0)]
# myscore_stab2 <- SelectionProportions(stab2)

######### run stab1 with ms = 1
stab1MS1 <- sharp::VariableSelection(
  xdata = simul$xdata, ydata = simul$ydata,
  implementation = CART1,
  family = "gaussian",
  maxsurrogate = 1)

myvars_stab1MS1 <- names(SelectedVariables(stab1MS1))[which(SelectedVariables(stab1MS1)!=0)]

####### run stab2 with ms = 1
stab2MS1 <- sharp::VariableSelection(
  xdata = simul$xdata, ydata = simul$ydata,
  implementation = CART2,
  family = "gaussian",
  Lambda = Lambda,
  maxsurrogate = 1)

myvars_stab2MS1 <- names(SelectedVariables(stab2MS1))[which(SelectedVariables(stab2MS1)!=0)]

######### run stab1 with ms = 2
stab1MS2 <- sharp::VariableSelection(
  xdata = simul$xdata, ydata = simul$ydata,
  implementation = CART1,
  family = "gaussian",
  maxsurrogate = 2)

myvars_stab1MS2 <- names(SelectedVariables(stab1MS2))[which(SelectedVariables(stab1MS2)!=0)]

####### run stab2 with ms = 2
stab2MS2 <- sharp::VariableSelection(
  xdata = simul$xdata, ydata = simul$ydata,
  implementation = CART2,
  family = "gaussian",
  Lambda = Lambda,
  maxsurrogate = 2)

myvars_stab2MS2 <- names(SelectedVariables(stab2MS2))[which(SelectedVariables(stab2MS2)!=0)]

####### Stability-lasso
lasso <- sharp::VariableSelection(
  xdata = simul$xdata, ydata = simul$ydata,
  implementation = PenalisedRegression,
  family = "gaussian")

myvars_LASSO <- names(SelectedVariables(lasso))[which(SelectedVariables(lasso)!=0)]

########## metrics
Methods <- c("CART", "stab1", "stab2", "stab1MS1", "stab2MS1", "stab1MS2", "stab2MS2", "LASSO")
Metrics <- vector("list", length = 8)
for (i in seq_along(Metrics)) {
  Metrics[[i]] <- getMetrics(getCM(Methods[i]))
}
Metrics <- as.data.frame(list.rbind(Metrics))
Metrics$Method <- Methods

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

p1 <- ggplot(df, aes(x=S, y=F1)) + geom_point() + ggtitle("stab")

################################ S vs F1 #######################################
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

p2 <- ggplot(df, aes(x=S, y=F1)) + geom_point() + ggtitle("stab2")

########### S vs F1 for stab3
F1 <- matrix(0, nrow=nrow(stab3$S_2d), ncol=ncol(stab3$S_2d))
for (i in 1:nrow(stab3$S_2d)) {
  for (j in 1:ncol(stab3$S_2d)) {
    myvars_stab3 <- Stable(stab3, argmax_id = c(i,j))
    myvars_stab3 <- names(myvars_stab3)[myvars_stab3 == 1]
    f1_stab3 <- getMetrics(getCM("stab3"))[3]
    F1[i,j] <- f1_stab3
  }
}

F1 <- c(F1)
S <- c(stab3$S_2d)
df <- data.frame(S = S, F1 = F1)
df <- df[!is.nan(df$S),]

p3 <- ggplot(df, aes(x=S, y=F1)) + geom_point() + ggtitle("stab-ms1")

########### S vs F1 for stab4
F1 <- matrix(0, nrow=nrow(stab4$S_2d), ncol=ncol(stab4$S_2d))
for (i in 1:nrow(stab4$S_2d)) {
  for (j in 1:ncol(stab4$S_2d)) {
    myvars_stab4 <- Stable(stab4, argmax_id = c(i,j))
    myvars_stab4 <- names(myvars_stab4)[myvars_stab4 == 1]
    f1_stab4 <- getMetrics(getCM("stab4"))[3]
    F1[i,j] <- f1_stab4
  }
}

F1 <- c(F1)
S <- c(stab4$S_2d)
df <- data.frame(S = S, F1 = F1)
df <- df[!is.nan(df$S),]

p4 <- ggplot(df, aes(x=S, y=F1)) + geom_point() + ggtitle("stab2-ms2")

########### S vs F1 for stab5
F1 <- matrix(0, nrow=nrow(stab5$S_2d), ncol=ncol(stab5$S_2d))
for (i in 1:nrow(stab5$S_2d)) {
  for (j in 1:ncol(stab5$S_2d)) {
    myvars_stab5 <- Stable(stab5, argmax_id = c(i,j))
    myvars_stab5 <- names(myvars_stab5)[myvars_stab5 == 1]
    f1_stab5 <- getMetrics(getCM("stab5"))[3]
    F1[i,j] <- f1_stab5
  }
}

F1 <- c(F1)
S <- c(stab5$S_2d)
df <- data.frame(S = S, F1 = F1)
df <- df[!is.nan(df$S),]

p5 <- ggplot(df, aes(x=S, y=F1)) + geom_point() + ggtitle("stab-ms2")

########### S vs F1 for stab6
F1 <- matrix(0, nrow=nrow(stab6$S_2d), ncol=ncol(stab6$S_2d))
for (i in 1:nrow(stab6$S_2d)) {
  for (j in 1:ncol(stab6$S_2d)) {
    myvars_stab6 <- Stable(stab6, argmax_id = c(i,j))
    myvars_stab6 <- names(myvars_stab6)[myvars_stab6 == 1]
    f1_stab6 <- getMetrics(getCM("stab6"))[3]
    F1[i,j] <- f1_stab6
  }
}

F1 <- c(F1)
S <- c(stab6$S_2d)
df <- data.frame(S = S, F1 = F1)
df <- df[!is.nan(df$S),]

p6 <- ggplot(df, aes(x=S, y=F1)) + geom_point() + ggtitle("stab2-ms2")

######## plots
ggarrange(p1,p3,p5,p2,p4,p6)
