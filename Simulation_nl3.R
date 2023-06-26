#####################
library(rpart)
library(sharp)
library(rlist)
library(parallel)
library(pbapply)
source("Functions.R")
source("stab1.R")
source("stab2.R")
source("SimulateX.R")
#####################
Simulation_study <- function(seed, n, pk, association, sd) {
  # set seed
  set.seed(seed)
  # simulate data
  X <- matrix(rnorm(n*pk),n,pk)
  mu1 <- X[,1]*X[,2]
  Y <- SimulateAssoc(mu1, association) + rnorm(n,0,sd)
  colnames(X) <- paste0("var", 1:pk)
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
  
  ####### run Stability-lasso
  lasso <- sharp::VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    implementation = PenalisedRegression,
    family = "gaussian")
  
  myvars <- names(SelectedVariables(lasso))[which(SelectedVariables(lasso)!=0)]
  selvars$LASSO <- myvars
  metrics$LASSO <- getMetrics(getCM(vars = myvars, Tvars = myvars_T, Fvars = myvars_F))
  
  ####### Cleanup output
  metrics <- as.data.frame(list.rbind(metrics))
  metrics$model <- rownames(metrics)
  
  ####### output
  return(list(metrics = metrics, selvars = selvars))
}

################ Try
# a <- Simulation_study(seed = 1, n = 1000, pk = 500, association = "Sine", sd = 1)
# b <- Simulation_study(seed = 1, n = 1000, pk = 500, association = "Linear", sd = 1)

################
# out <- vector("list", length=length(myassoc))
# for (i in 1:length(myassoc)) {
#   out[[i]] <- Simulation_study(seed = 1, association = myassoc[i], 
#                    n = n, pk = pk, sd = sd)
# }


################ test all associations
n <- 1000
pk <- 500
sd <- 1
nchunks <- 5
myassoc <- c("Heteroscedastic", "Absolute", "Quadratic", "Cosine", "Sine", "Indicator",
             "SquareWave", "Sawtooth", "Exponential", "Sigmoidal", "Cubic", "Linear"  )

no_cores <- nchunks
cl <- makeCluster(no_cores)

clusterEvalQ(cl, library(rpart))
clusterEvalQ(cl, library(sharp))
clusterEvalQ(cl, library(fake))
clusterEvalQ(cl, library(rlist))
clusterExport(cl, c("myassoc", "n", "pk", "sd",
                    "CART1", "CART2", "getBeta", "getCM", "getMetrics", "Simulation_study",
                    "SimulateAssoc"))

out <- pblapply(1:length(myassoc),
                function(i) {Simulation_study(seed = 1, association = myassoc[i], 
                                              n = n, pk = pk, sd = sd)},
                cl = cl)

stopCluster(cl)

############ Cleanup output
metrics <- vector("list", length=length(myassoc))
for (i in seq_along(out)) {
  metrics[[i]] <- out[[i]]$metrics
}
names(metrics) <- myassoc
metrics <- as.data.frame(list.rbind(metrics))
metrics$Relationship <- rep(myassoc, each=4)

selvars <- vector("list", length=length(myassoc))
for (i in seq_along(out)) {
  selvars[[i]] <- out[[i]]$selvars
}
names(selvars) <- myassoc

############ Visualise
library(ggplot2)
library(forcats)
metrics$model <- factor(metrics$model, levels = unique(metrics$model))
metrics$Relationship <- factor(metrics$Relationship, levels=myassoc)
title <- paste0("Y = nl(x1*x2), numrep = 1, n = 1000, pk = 500, sd = ", sd)

p1 <- ggplot(metrics, aes(y=Recall, x=Relationship, fill=forcats::fct_rev(model))) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "PuOr", 
                     guide = guide_legend(reverse = TRUE),
                     name = element_blank()) +
  coord_flip() +
  theme_grey() +
  theme(axis.title.y = element_blank())

p2 <- ggplot(metrics, aes(y=Precision, x=Relationship, fill=forcats::fct_rev(model))) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "PuOr", 
                     guide = guide_legend(reverse = TRUE),
                     name = element_blank()) +
  coord_flip() +
  theme_grey() +
  theme(axis.title.y = element_blank())

p3 <- ggplot(metrics, aes(y=F1, x=Relationship, fill=forcats::fct_rev(model))) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "PuOr", 
                     guide = guide_legend(reverse = TRUE),
                     name = element_blank()) +
  coord_flip() +
  theme_grey() +
  theme(axis.title.y = element_blank())

p4 <- ggplot(metrics, aes(y=Specificity, x=Relationship, fill=forcats::fct_rev(model))) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "PuOr", 
                     guide = guide_legend(reverse = TRUE),
                     name = element_blank()) +
  coord_flip() +
  theme_grey() +
  theme(axis.title.y = element_blank())

library(ggpubr)
p5 <- ggarrange(p1,p2,p3,p4, 
          common.legend = TRUE,
          legend = "bottom")

annotate_figure(p5, top = text_grob(title, color = "black", face = "bold", size = 14))

##### save output
x1_mult_x2_sd1 <- list(metrics = metrics, selvars = selvars)
saveRDS(x1_mult_x2_sd1, "x1_plus_x2_sd1.rds")
