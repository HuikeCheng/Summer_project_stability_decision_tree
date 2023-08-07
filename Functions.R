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
  if (is.nan(precision)) {precision <- 0}
  f1 <- 2*recall*precision/(recall+precision)
  if (is.nan(f1)) {f1 <- 0}
  mcc <- (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  # output
  metrics <- c(hd, f1, recall, precision, specificity, mcc)
  names(metrics) <- c("HD", "F1", "Recall", "Precision", "Specificity", "MCC")
  return(metrics)
}

selvars <- c(0,0,1,1,0,0)
theta <- c(0,1,0,1,0,0)
getMetrics(selvars, theta)

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

