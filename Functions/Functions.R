# Functions

getMetrics <- function(selvars, theta) {
  # Hamming distance
  hd <- as.numeric(sum(selvars != theta))
  # confusion matrix
  TP <- as.numeric(sum((selvars+theta) == 2))
  TN <- as.numeric(sum((selvars+theta) == 0))
  FP <- as.numeric(sum((selvars-theta) == 1))
  FN <- as.numeric(sum((selvars-theta) == -1))
  # other metrics
  specificity <- TN/(TN+FP)
  recall <- TP/(TP+FN)
  precision <- TP/(TP+FP)
  if (is.nan(precision)) {
    f1 <- NaN
  } else {
    f1 <- 2*recall*precision/(recall+precision)
  }
  mcc <- (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  # output
  metrics <- c(hd, f1, recall, precision, specificity, mcc)
  names(metrics) <- c("HD", "F1", "Recall", "Precision", "Specificity", "MCC")
  return(metrics)
}

# selvars <- c(0,0,1,1,0,0)
# theta <- c(0,1,0,1,0,0)
# getMetrics(selvars, theta)

HugeAdjacency <- function(pk = 10, topology = "random", nu = 0.1, ...) {
  # Storing extra arguments
  extra_args <- list(...)
  
  # Extracting relevant extra arguments
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = huge::huge.generator)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("n", "d", "prob", "graph", "verbose")]
  
  # Running simulation model
  mymodel <- do.call(huge::huge.generator, args = c(
    list(
      n = 2, d = sum(pk), prob = nu,
      graph = topology, verbose = FALSE
    ),
    tmp_extra_args
  ))
  theta <- as.matrix(mymodel$theta)
  
  # Re-organising the variables to avoid having centrality related to variable ID (e.g. for scale-free models)
  ids <- sample(ncol(theta))
  theta <- theta[ids, ids]
  
  return(theta)
}

############# ZScore over varimp
ZScore <- function(varimp) {
  vithr <- seq(min(varimp), max(varimp), length.out = 100)
  out <- vector("numeric", length = length(vithr))
  for (i in 1:length(vithr)) {
    G1 <- varimp[which(varimp <= vithr[i])]
    G2 <- varimp[which(varimp > vithr[i])]
    
    out[i] <- (mean(G2)-mean(G1))/sqrt((var(G2)/length(G2))+(var(G1)/length(G1)))
    
  }
  return(vithr[which.max(out)])
}


# varimp <- rf$importance
# output <- ZScore(varimp)

#######################
getImpRF_new<-function(x,y,ntree=500,num.trees=ntree,...){
  if(inherits(y,"Surv")){
    x$shadow.Boruta.time<-y[,"time"]
    x$shadow.Boruta.status<-y[,"status"]
    return(ranger::ranger(data=x,
                          dependent.variable.name="shadow.Boruta.time",
                          status.variable.name="shadow.Boruta.status",
                          num.trees=num.trees,importance="permutation",
                          scale.permutation.importance=TRUE,
                          write.forest=FALSE,
                          num.threads = 1,...)$variable.importance)
  }
  #Abusing the fact that Boruta disallows attributes with names
  # starting from "shadow"
  x$shadow.Boruta.decision<-y
  ranger::ranger(data=x,dependent.variable.name="shadow.Boruta.decision",
                 num.trees=num.trees,importance="permutation",
                 scale.permutation.importance=TRUE,
                 write.forest=FALSE,
                 num.threads = 1,...)$variable.importance
}
#comment(getImpRF_new)<-'ranger normalized permutation importance'

# getCM <- function(vars, Tvars, Fvars) {
#   cm <- matrix(0, nrow = 2, ncol = 2)
#   cm[1,1] <- sum(Tvars %in% vars)
#   cm[1,2] <- sum(!(Tvars %in% vars))
#   cm[2,1] <- sum(Fvars %in% vars)
#   cm[2,2] <- sum(!(Fvars %in% vars))
#   return(cm)
# }

# Tvars <- c("x1", "x2", "x3")
# Fvars <- c("x4", "x5", "x6")
# vars <- c("x1", "x2")
# getCM(vars, Tvars, Fvars)

# getMetrics <- function(cm) {
#   recall <- cm[1,1]/(cm[1,1] + cm[1,2])
#   precision <- cm[1,1]/(cm[1,1] + cm[2,1])
#   f1 <- 2*recall*precision/(recall + precision)
#   specificity <- cm[2,2]/(cm[2,2] + cm[2,1])
#   
#   if (is.nan(precision)) {precision <- 0}
#   if (is.nan(f1)) {f1 <- 0}
#   
#   metrics <- c(recall, precision, f1, specificity)
#   names(metrics) <- c("Recall", "Precision", "F1", "Specificity")
#   return(metrics)
# }

# cm <- getCM("cart")
# getMetrics(cm)

# getBeta <- function(vars, Tvars, beta) {
#   TP <- Tvars[Tvars %in% vars]
#   FN <- Tvars[!(Tvars %in% vars)]
#   meanTP <- mean(abs(beta[rownames(beta) %in% TP]))
#   meanFN <- mean(abs(beta[rownames(beta) %in% FN]))
#   minTP <- min(abs(beta[rownames(beta) %in% TP]))
#   minFN <- min(abs(beta[rownames(beta) %in% FN]))
#   maxTP <- max(abs(beta[rownames(beta) %in% TP]))
#   maxFN <- max(abs(beta[rownames(beta) %in% FN]))
#   out <- c(minTP, meanTP, maxTP, minFN, meanFN, maxFN)
#   names(out) <- c("minTP", "meanTP", "maxTP", "minFN", "meanFN", "maxFN")
#   if (any(is.nan(out))) {
#     print("No FN, set NaN and Inf to 0")
#     out[is.nan(out)] <- 0
#     out[is.infinite(out)] <- 0
#   }
#   return(out)
# }

Bin2Dec <- function(x) {
  x <- sum(2^(which(rev(x) == 1) - 1))
  return(x)
}

## v1
# reorderCV <- function(varname, vp) {
#   order <- which(colSums(vp) != 0) # which levels contain the variable
#   order <- order[length(order)] # the the last level with the variable
#   for(i in order:1) {
#     if (sum(vp[,i] != 0) > 0) {
#       row_num <- which(vp[,i] != 0)[1]
#       previous <- which(vp[row_num, 1:i-1] != 0)
#       if (length(previous) > 0) {
#         previous <- previous[length(previous)]
#         previous_path <- vp[row_num, previous]
#         if (previous_path == 1) {
#           order <- c(order, previous)
#         } else if (previous_path == 2) {
#           order <- c(previous, order)
#         }
#       } else {
#         order <- c(order, i)
#       }
#     }
#   }
#   order <- unique(order)
#   return(order)
# }

### v2
# reorderCV <- function(vp) {
#   layers <- which(colSums(vp) != 0) # which layers contain the variable
#   layers <- rev(layers)
#   order <- layers[1]
#   for(i in seq_along(layers)) {
#     current <- layers[i]
#     row_num <- which(vp[,layers[i]] != 0)[1] # get the row number of this variable in vp, 
#     # only take the 1st index, as in current SimulateTree, vars can only occur once in each layer, 
#     # so all entries in a layer has same ancestors in previous layers
#     previous <- layers[-c(1:i)]
#     previous <- previous[vp[row_num, previous] != 0]
#     if (length(previous) > 0) {
#       previous <- previous[1]
#       previous_path <- vp[row_num, previous]
#       if (!current %in% order) {
#         if (previous %in% order) {
#           if (previous_path == 1) {
#             order <- c(current, order)
#           } else if (previous_path == 2) {
#             order <- c(order, current)
#           }
#         } else {
#           order <- c(order, current)
#         }
#       } else {
#         if (previous_path == 1) {
#           order <- c(order, previous)
#         } else if (previous_path == 2) {
#           order <- c(previous, order)
#         }
#       }
#     } else if (i+1 <= length(layers)){
#         order <- c(order, layers[i+1])
#       }
#   #order <- unique(order)
#   }
#   return(order)
# }

#reorderCV(vp)

# still has problems: consider cases where share same ancestor at root but all just have this one ancestor
# sometimes the current will be missed
# consider case where one has ancestor so order skipped the middle, and the middle one does not have ancestor,
# then the middle would be missed

## v3
reorderCV <- function(vp) {
  layers <- which(colSums(vp) != 0) # get layers that contain the variable
  layers <- rev(layers) # start from last
  order <- layers[1]
  for(i in seq_along(layers)) {
    current <- layers[i]
    row_num <- which(vp[,layers[i]] != 0)[1]
    previous <- layers[-c(1:i)]
    previous <- previous[vp[row_num, previous] != 0]
    if (length(previous) > 0) {
      previous <- previous[1]
      previous_path <- vp[row_num, previous]
      if (!previous %in% order) {
        if(!current %in% order) {
          order <- c(order, current)
        }
        if (previous_path == 1) {
          order <- c(order, previous)
        } else {
          order <- c(previous, order)
        }
      } else {
        if (previous_path == 1) {
          order <- c(current, order)
        } else {
          order <- c(order, current)
        }
      }
    } else {
      if(!current %in% order) {
        order <- c(order, current)
      }
    }
  }
  order <- unique(order)
  return(order)
}

#######################
rfvimptest_new <- function(data, yname, Mmax = 500, varnames = NULL, p0 = 0.06, p1 = 0.04, alpha = 0.05, beta = 0.2, A = 0.1, B = 10, h = 8, nperm = 1,
                           ntree = 500,
                           progressbar = TRUE,
                           test = c("general", "twosample")[1],
                           type = c("SPRT", "SAPT", "pval", "certain", "complete")[1], condinf=FALSE, ...) {
  starttime <- Sys.time()
  # @seealso \code{\link{predict.divfor}}
  
  if(any(is.na(data))) {
    missvariables <- paste(names(data)[apply(data, 2, function(x) any(is.na(x)))], collapse = ", ")
    stop(paste0("Missing data in columns: ", missvariables, ". Please provide complete data or consider setting condinf=TRUE."))
  }
  
  if(!condinf & test == "twosample")
    stop("'twosample' approach only useable if condinf=TRUE.")
  
  if(!condinf & nperm > 1)
    stop("Values of 'nperm' different than 1 are only useable if condinf=TRUE.")
  
  if (progressbar) pb <- txtProgressBar(min = 1, max = Mmax, initial = 1, width = 10, style = 3, char = "|")
  
  if (type == "SPRT") {
    A <- beta / (1 - alpha)
    B <- (1 - beta) / alpha
  }
  if (type %in% c("SPRT", "SAPT")) {
    logA <- log(A)
    logB <- log(B)
    help1 <- log((1 - p0) / (1 - p1))
    help2 <- log((p1 * (1 - p0)) / (p0 * (1 - p1)))
  }
  
  stop_crits <- switch(type,
                       SPRT = list((logA + 1:Mmax * help1) / help2, (logB + 1:Mmax * help1) / help2),
                       SAPT = list((logA + 1:Mmax * help1) / help2, (logB + 1:Mmax * help1) / help2),
                       pval = list(rep(h, times = Mmax), rep(h, times = Mmax)),
                       certain = list(rep(alpha*Mmax, times = Mmax), Mmax*alpha - Mmax + 1:Mmax),
                       complete = NULL)
  
  if (is.null(varnames))
    varnames <- names(data)[names(data)!=yname]
  
  if (!condinf) {
    rfmod <- ranger::ranger(data = data, dependent.variable.name=yname, num.tree=ntree, importance = "permutation", num.threads = 1)
    vimp_orig <- list()
    vimp_orig$values <- rfmod$variable.importance[varnames]
  }
  else {
    rfmod <- party::cforest(as.formula(paste(yname, " ~ .", sep="")), data = data, controls = party::cforest_unbiased(ntree=ntree, ...))
    vimp_orig <- permimp::permimp(rfmod, whichxnames = varnames, nperm = nperm, asParty = TRUE, progressBar = FALSE)
  }
  
  
  if (test == "general") {
    testresult <-
      lapply(varnames, function(v) {
        permdata <- data
        permvimps <- c()
        for (m in 1:Mmax) {
          if (progressbar) {setTxtProgressBar(pb, m)}
          if (m == 1 & progressbar) cat(" of variable", v)
          permdata[, v] <- sample(permdata[, v])
          
          if (!condinf) {
            permmod <- ranger::ranger(data = permdata, dependent.variable.name=yname, num.tree=ntree, importance="permutation", num.threads = 1)
            permvimps <- c(permvimps, permmod$variable.importance[v])
          }
          else {
            permmod <- party::cforest(as.formula(paste(yname, " ~ .", sep="")), data = permdata, controls = party::cforest_unbiased(ntree=ntree, ...))
            permvimps <- c(permvimps, permimp::permimp(permmod, whichxnames = v, nperm = nperm, asParty = TRUE, progressBar = FALSE)$values)
          }
          d <- sum(permvimps >= vimp_orig$values[v])
          if (type == "certain") {
            if (d > stop_crits[[1]][m]) {result <- "keep H0"; pvalue <- NA; break} else if (d <= stop_crits[[2]][m]) {result <- "accept H1"; pvalue <- NA; break}
          } else if (type %in% c("SPRT", "SAPT")) {
            if (d >= stop_crits[[1]][m]) {result <- "keep H0"; pvalue <- NA; break} else if (d <= stop_crits[[2]][m]) {result <- "accept H1"; pvalue <- NA; break}
          } else if (type == "pval") {
            if (d == h) {
              pvalue <- d/m
              result <- ifelse(pvalue > 0.05, "keep H0", "accept H1")
              break
            }
          }
        }
        if (m == Mmax) {
          if (type == "pval") {
            if (d < h) {
              pvalue <- (d + 1) / (Mmax + 1)
              result <- ifelse(pvalue > 0.05, "keep H0", "accept H1")
            }
            else {
              pvalue <- d / Mmax
              result <- ifelse(pvalue > 0.05, "keep H0", "accept H1")
            }
          } else  if (type == "complete") {
            pvalue <- d / Mmax
            result <- ifelse(pvalue > 0.05, "keep H0", "accept H1")
          } else {
            pvalue <- NA
            result <- ifelse(d / Mmax > 0.05, "keep H0", "accept H1")
          }
        }
        if (progressbar) cat(" - finished", "\n")
        list(testres=result, pvalue=pvalue, stoppedearly=ifelse(m < Mmax, "yes", "no"), permperf=m)
      })
  } else if (test == "twosample") {
    testresult <-
      lapply(varnames, function(v) {
        permdata <- data
        permdata[, v] <- sample(permdata[, v])
        permmod <- party::cforest(as.formula(paste(yname, " ~ .", sep="")), data = permdata, controls = party::cforest_unbiased(ntree=ntree, ...))
        permmodvimps <- permimp::permimp(permmod, whichxnames = v, nperm = nperm, asParty = TRUE, progressBar = FALSE)
        permvimps <- c()
        for (m in 1:Mmax) {
          if (progressbar) {setTxtProgressBar(pb, m)}
          if (m == 1 & progressbar) cat(" of variable", v)
          permvimps <- c(permvimps, mean(c(vimp_orig$perTree[, v], permmodvimps$perTree[, v])[sample(1:(2*ntree), ntree)]))
          d <- sum(permvimps >= vimp_orig$values[v])
          if (type == "certain") {
            if (d > stop_crits[[1]][m]) {result <- "keep H0"; pvalue <- NA; break} else if (d <= stop_crits[[2]][m]) {result <- "accept H1"; pvalue <- NA; break}
          } else if (type %in% c("SPRT", "SAPT")) {
            if (d >= stop_crits[[1]][m]) {result <- "keep H0"; pvalue <- NA; break} else if (d <= stop_crits[[2]][m]) {result <- "accept H1"; pvalue <- NA; break}
          } else if (type == "pval") {
            if (d == h) {
              pvalue <- d/m
              result <- ifelse(pvalue > 0.05, "keep H0", "accept H1")
              break
            }
          }
        }
        if (m == Mmax) {
          if (type == "pval") {
            if (d < h) {
              pvalue <- (d + 1) / (Mmax + 1)
              result <- ifelse(pvalue > 0.05, "keep H0", "accept H1")
            }
            else {
              pvalue <- d / Mmax
              result <- ifelse(pvalue > 0.05, "keep H0", "accept H1")
            }
          } else  if (type == "complete") {
            pvalue <- d / Mmax
            result <- ifelse(pvalue > 0.05, "keep H0", "accept H1")
          } else {
            pvalue <- NA
            result <- ifelse(d / Mmax > 0.05, "keep H0", "accept H1")
          }
        }
        if (progressbar) cat(" - finished", "\n")
        list(testres=result, pvalue=pvalue, stoppedearly=ifelse(m < Mmax, "yes", "no"), permperf=m)
      })
  }
  stoptime <- Sys.time()
  time_elapsed <- paste0(round(as.numeric(difftime(stoptime, starttime, units = "secs")), 1), " seconds")
  testres <- sapply(testresult, function(x) x$testres)
  pvalues <- sapply(testresult, function(x) x$pvalue)
  stoppedearly <- sapply(testresult, function(x) x$stoppedearly)
  perms <- sapply(testresult, function(x) x$permperf)
  names(testres) <- names(pvalues) <- names(stoppedearly) <- names(perms) <- names(vimp_orig$values)
  
  result <- list(testtype = paste(test, type, sep=", "), varimp=vimp_orig$values, testres = testres, pvalues = pvalues,
                 stoppedearly = stoppedearly, perms = perms, Mmax=Mmax, ntree=ntree, comptime = time_elapsed)
  
  class(result) <- "rfvimptest"
  
  return(result)
  
}

#######################
getImpXgboost_new <-function(x,y,nrounds=5,verbose=0,nthread=1){
  for(e in 1:ncol(x)) x[,e]<-as.numeric(x[,e])
  xgboost::xgb.importance(
    model=xgboost::xgboost(
      data=as.matrix(x),
      label=y,
      nrounds=nrounds,
      verbose=verbose,
      nthread=nthread
    )
  )->imp
  stats::setNames(rep(0,ncol(x)),colnames(x))->ans
  ans[imp$Feature]<-imp$Gain
  ans
}
