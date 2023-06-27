#####################
library(rpart)
library(sharp)
library(rlist)
source("Functions.R")
source("stab1.R")
source("stab2.R")
source("SimulateX.R")
#####################
n <- 1000
pk <- 500
X <- matrix(rnorm(n*pk),n,pk)
R <- 1:2
beta <- rep(1, 2)
mu1 <- X[,R]%*%beta
Y <- SimulateAssoc(mu1, "Sawtooth") + rnorm(n,0,1)
colnames(X) <- paste0("var", 1:500)
simul <- list(ydata=Y, xdata=X)

# vars
myvars_T <- paste0("var", c(1,2))
myvars_F <- paste0("var", c(3:500))

# set up output
metrics <- vector("list", length = 4)
selvars <- vector("list", length = 4)
names(metrics) <- c("CART", "stab1", "stab2", "LASSO")
names(selvars) <- c("CART", "stab1", "stab2", "LASSO")

######## run CART
mydata <- data.frame(ydata = simul$ydata[, 1], simul$xdata)
mycart <- rpart(ydata ~ ., mydata, method="anova")

myvars <- unique(rownames(mycart$splits))
selvars$CART <- myvars
metrics$CART <- getMetrics(getCM(vars = myvars, Tvars = myvars_T, Fvars = myvars_F))

######### run stab1
stab1 <- sharp::VariableSelection(
  xdata = simul$xdata, ydata = simul$ydata,
  implementation = CART1,
  family = "gaussian",
  maxsurrogate = 0)

myvars <- names(SelectedVariables(stab1))[which(SelectedVariables(stab1)!=0)]
selvars$stab1 <- myvars
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
selvars$stab2 <- myvars
metrics$stab2 <- getMetrics(getCM(vars = myvars, Tvars = myvars_T, Fvars = myvars_F))

####### Stability-lasso
lasso <- sharp::VariableSelection(
  xdata = simul$xdata, ydata = simul$ydata,
  implementation = PenalisedRegression,
  family = "gaussian")

myvars <- names(SelectedVariables(lasso))[which(SelectedVariables(lasso)!=0)]
selvars$LASSO <- myvars
metrics$LASSO <- getMetrics(getCM(vars = myvars, Tvars = myvars_T, Fvars = myvars_F))

############ Cleanup output
metrics <- as.data.frame(list.rbind(metrics))
metrics$model <- rownames(metrics)


