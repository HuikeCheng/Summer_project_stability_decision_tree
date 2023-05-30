######### load packages
library(rpart)
library(sharp)
library(fake)

######### simulate
set.seed(1)
simul <- SimulateRegression(pk = 50)

######### true vars
myvars_T <- rownames(simul$beta)[which(simul$beta != 0)]
myvars_F <- rownames(simul$beta)[which(simul$beta == 0)]

######### run stab
stab <- sharp::VariableSelection(
  xdata = simul$xdata, ydata = simul$ydata,
  implementation = CART,
  family = "gaussian",
)
myvars_stab <- names(SelectedVariables(stab))[which(SelectedVariables(stab)!=0)]
myscore_stab <- SelectionProportions(stab)

######## run CART
mydata <- data.frame(ydata = simul$ydata[, 1], simul$xdata)
mycart <- rpart(ydata ~ ., mydata, method="anova", )
myvars_cart <- names(mycart$variable.importance)
myscore_cart <- mycart$variable.importance

####### Confusion matrix
TP_stab <- sum(myvars_T %in% myvars_stab)
FP_stab <- sum(myvars_F %in% myvars_stab)
TN_stab <- sum(!(myvars_F %in% myvars_stab))
FN_stab <- sum(!(myvars_T %in% myvars_stab))

TP_cart <- sum(myvars_T %in% myvars_cart)
FP_cart <- sum(myvars_F %in% myvars_cart)
TN_cart <- sum(!(myvars_F %in% myvars_cart))
FN_cart <- sum(!(myvars_T %in% myvars_cart))

########## metrics
sensitivity_stab <- TP_stab/length(myvars_T)
precision_stab <- TP_stab/length(myvars_stab)
f1_stab <- 2*sensitivity_stab*precision_stab/(sensitivity_stab+precision_stab)

sensitivity_cart <- TP_cart/length(myvars_T)
precision_cart <- TP_cart/length(myvars_cart)
f1_cart <- 2*sensitivity_cart*precision_cart/(sensitivity_cart+precision_cart)
