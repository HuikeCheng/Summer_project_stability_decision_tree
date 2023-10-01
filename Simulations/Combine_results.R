# load packages
library(rlist)
library(tidyverse)
library(parallel)
library(pbapply)

# parameters involved
params <- read.delim("params_c2.txt", sep = " ", header = FALSE)
params <- params[,-c(1, ncol(params))]

scenario <- "Complex2"

# read in datasets
out <- vector("list", length=nrow(params))
for (i in seq_along(out)) {
  folder <- paste0(c(scenario, params[i,]), collapse="_")
  out[[i]] <- readRDS(paste("sCART1Outputs",folder,"Metrics.rds", sep="/"))
}

# combine
out <- list.rbind(out)
out$relationship <- rep(params[,1], each = 40000)
out$n <- rep(params[,2], each = 40000)
out$pk <- rep(params[,3], each = 40000)
out$ev_xy <- rep(params[,4], each = 40000)
out$numpair <- rep(params[,5], each = 40000)
out$numvar <- rep(params[,6], each = 40000)

# save
#saveRDS(out, "Mult.rds")

############## Tree

# parameters involved
params <- read.delim("params_c2.txt", sep = " ", header = FALSE)
params <- params[,-c(1, 2, ncol(params))]

scenario <- "Complex2"

# read in datasets
out <- vector("list", length=nrow(params))
for (i in seq_along(out)) {
  folder <- paste0(c(scenario, params[i,]), collapse="_")
  out[[i]] <- readRDS(paste("sCART1Outputs",folder,"Metrics.rds", sep="/"))
}

# combine
out <- list.rbind(out)
out$n <- rep(params[,1], each = 1000)
out$pk <- rep(params[,2], each = 1000)
out$ev_xy <- rep(params[,3], each = 1000)

# save
#saveRDS(out, "C2_sCART1_1000.rds")


#################### Sels #####################
selWS <- function(theta, mysel) {
  tmp <- theta + mysel
  if (sum(tmp==2) == sum(theta)) {
    return(1)
  } else {
    return(0)
  }
}
# parameters involved
params <- read.delim("params_c2.txt", sep = " ", header = FALSE)
params <- params[,-c(1, ncol(params))]

scenario <- "Complex2"

# read in datasets
out <- vector("list", length=nrow(params))
for (i in seq_along(out)) {
  folder <- paste0(c(scenario, params[i,]), collapse="_")
  df <- readRDS(paste("sCART1Outputs",folder,"Sels.rds", sep="/"))
  theta <- readRDS(paste("sCART1Outputs",folder,"Theta.rds", sep="/"))
  numsel <- rowSums(df[,-ncol(df)])
  numWS <- apply(df[,-ncol(df)], 1, function(x) selWS(theta, x))
  sel <- data.frame(numsel = numsel, WS = numWS, model = df[,ncol(df)])
  out[[i]] <- sel
}

# combine
out <- list.rbind(out)
out$relationship <- rep(params[,1], each = 40000)
out$n <- rep(params[,2], each = 40000)
out$pk <- rep(params[,3], each = 40000)
out$ev_xy <- rep(params[,4], each = 40000)
out$numpair <- rep(params[,5], each = 40000)
out$numvar <- rep(params[,6], each = 40000)

# save
#saveRDS(out, "C2Sel_sCART1_1000.rds")

#############################################################################
################################# Stability #################################
#############################################################################
getJaccard <- function(set1, set2) {
  tmp <- set1 + set2
  intersect <- sum(tmp==2)
  union <- sum(tmp > 0)
  return(intersect/union)
}

getStab <- function(df) {
  numrep <- ncol(df)
  stab <- matrix(NA, nrow=numrep, ncol=numrep)
  for (i in 1:(numrep-1)) {
    stab[(i+1):numrep,i] <- apply(data.frame((i+1):numrep), 1, 
                                  function(x) getJaccard(df[,i], df[,x]))
  }
  return(mean(stab, na.rm = TRUE))
}

runStab <- function(scenario, params, models) {
  folder <- paste0(c(scenario, params), collapse="_")
  df <- readRDS(paste("Outputs",folder,"Sels.rds", sep="/"))
  numModel <- length(models)
  tmp <- vector("numeric", length=numModel)
  
  for (i in seq_along(models)) {
    df_md <- df %>% filter(model==models[i]) %>% select(-model)
    df_md <- t(df_md)
    tmp[i] <- getStab(df_md) 
  }
  
  return(data.frame(stab = tmp, model = models,
                    n = rep(params[1], numModel),
                    pk = rep(params[2], numModel),
                    ev_xy = rep(params[3], numModel)))
}

# df in getStab should have simulations as columns and vars as rows
# df_sCART1 <- df %>% filter(model=="sCART1") %>% select(-model)
# df_sCART1 <- t(df_sCART1)

########################## Complex 2 #######################################
params <- read.delim("params_c2.txt", sep = " ", header = FALSE)
params <- params[,-c(1, ncol(params))]

scenario <- "Complex2"
models <- c("CART", "sCART2", "sLASSO")

# read in datasets
out <- vector("list", length=nrow(params))
for (i in 1:nrow(params)) {
  folder <- paste0(c(scenario, params[i,]), collapse="_")
  df <- readRDS(paste("Outputs",folder,"Sels.rds", sep="/"))
  numModel <- length(models)
  tmp <- vector("numeric", length=numModel)
  for (j in seq_along(models)) {
    df_md <- df %>% filter(model==models[j]) %>% select(-model)
    df_md <- t(df_md)
    tmp[j] <- getStab(df_md)
  }
  out[[i]] <- data.frame(stab = tmp, model = models,
                         n = rep(params[i,1], numModel),
                         pk = rep(params[i,2], numModel),
                         ev_xy = rep(params[i,3], numModel))
}

# combine
out <- list.rbind(out)

#saveRDS(out, "stab_sCART1_C2_test.rds")

########### Parallelise
nchunks <- 5
no_cores <- nchunks
cl <- makeCluster(no_cores)

clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(rlist))
clusterExport(cl, c("scenario","params", "models","getJaccard", "getStab", "runStab"))

out <- pblapply(1:nrow(params),
                function(i) {runStab(scenario = scenario, params = params[i,], 
                                     models = models)},
                cl = cl)

stopCluster(cl)



####################### Test #######################
# stab <- matrix(NA, nrow=ncol(df_sCART1), ncol=ncol(df_sCART1))
# numsim <- ncol(stab)
# 
# start.time <- Sys.time()
# for (i in 1:(numsim-1)) {
#   stab[i,(i+1):numsim] <- apply(data.frame((i+1):numsim), 1, 
#                                 function(x) getJaccard(df_sCART1[,i], df_sCART1[,x]))
# }
# end.time <- Sys.time()
# time.taken_1 <- end.time-start.time #20.26279 s
# 
# start.time <- Sys.time()
# for (i in 1:(numsim-1)) {
#   stab[(i+1):numsim,i] <- apply(data.frame((i+1):numsim), 1, 
#                                 function(x) getJaccard(df_sCART1[,i], df_sCART1[,x]))
# }
# end.time <- Sys.time()
# time.taken_2 <- end.time-start.time # 20.15672 s
# 
# start.time <- Sys.time()
# for (i in 1:(numsim-1)) {
#   for (j in (i+1):numsim) {
#     stab[i,j] <- getJaccard(df_sCART1[,i], df_sCART1[,j])
#   }
# }
# end.time <- Sys.time()
# time.taken_3 <- end.time-start.time # 22.41977 s

################ update sCART1
summary(factor(Mult_sCART1_Stability$model))
TreeNR5Sel_10000 <- TreeNR5Sel_10000 %>% filter(model != "sCART1")
df <- rbind(TreeNR5Sel_10000, TreeNR5Sel_sCART1_10000)
#saveRDS(df, "../Continuous_10000_v2/TreeNR5Sel.rds")

colnames(C2_sCART1_Stability)[3:5] <- c("n", "pk", "ev_xy")
C2_Stability <- C2_Stability[,-c(4,5,7,8,10,11)]
colnames(C2_Stability)[3:5] <- c("n", "pk", "ev_xy")

Mult_Stab <- rbind(Mult_stability, Mult_sCART1_stability)
#saveRDS(Mult_Stab, "../Continuous_10000_v2/MultStab.rds")

MultSel_sCART1_10000 <- MultSel_sCART1_10000 %>% select(-numNoise)



# out <- rbind(Plus, Plus_4610)
# saveRDS(out, "../Continuous_10000_v3/Plus.rds")
# 
# out <- rbind(PlusSel, Plus4610Sel)
# saveRDS(out, "../Continuous_10000_v3/PlusSel.rds")
# 
# out <- rbind(PlusStab, Plus4610Stab)
# saveRDS(out, "../Continuous_10000_v3/PlusStab.rds")
