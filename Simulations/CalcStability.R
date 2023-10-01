############ packages
library(rlist)
library(tidyverse)
library(parallel)

############ functions
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
  df <- readRDS(paste("ReMult4Outputs",folder,"Sels.rds", sep="/"))
  numModel <- length(models)
  tmp <- vector("numeric", length=numModel)
  
  for (i in seq_along(models)) {
    df_md <- df %>% filter(model==models[i]) %>% select(-model)
    df_md <- t(df_md)
    tmp[i] <- getStab(df_md) 
  }
  
  return(data.frame(stab = tmp, model = models,
                    relationship = rep(params[,1], numModel),
                    n = rep(params[,2], numModel),
                    pk = rep(params[,3], numModel),
                    ev_xy = rep(params[,4], numModel),
                    numpair = rep(params[,5], numModel),
                    numvar = rep(params[,6], numModel)))
}

############ set up
params <- read.delim("params_mult_2.txt", sep = " ", header = FALSE)
params <- params[,-c(1, ncol(params))]

scenario <- "Mult"
models <- c("CART", "sCART1", "sCART2", "sLASSO")

########### Parallelise
no_cores <- 10
cl <- makeCluster(no_cores)

clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(rlist))
clusterExport(cl, c("scenario","params", "models","getJaccard", "getStab", "runStab"))

out <- parLapply(1:nrow(params),
                function(i) {runStab(scenario = scenario, params = params[i,], 
                                     models = models)},
                cl = cl)

stopCluster(cl)

############## Cleanup data
out <- list.rbind(out)

saveRDS(out, file = "Mult4_Stability.rds")