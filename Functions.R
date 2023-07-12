getCM <- function(vars, Tvars, Fvars) {
  cm <- matrix(0, nrow = 2, ncol = 2)
  cm[1,1] <- sum(Tvars %in% vars)
  cm[1,2] <- sum(!(Tvars %in% vars))
  cm[2,1] <- sum(Fvars %in% vars)
  cm[2,2] <- sum(!(Fvars %in% vars))
  return(cm)
}

# Tvars <- c("x1", "x2", "x3")
# Fvars <- c("x4", "x5", "x6")
# vars <- c("x1", "x2")
# getCM(vars, Tvars, Fvars)

getMetrics <- function(cm) {
  recall <- cm[1,1]/(cm[1,1] + cm[1,2])
  precision <- cm[1,1]/(cm[1,1] + cm[2,1])
  f1 <- 2*recall*precision/(recall + precision)
  specificity <- cm[2,2]/(cm[2,2] + cm[2,1])
  
  if (is.nan(precision)) {precision <- 0}
  if (is.nan(f1)) {f1 <- 0}
  
  metrics <- c(recall, precision, f1, specificity)
  names(metrics) <- c("Recall", "Precision", "F1", "Specificity")
  return(metrics)
}

# cm <- getCM("cart")
# getMetrics(cm)

getBeta <- function(vars, Tvars, beta) {
  TP <- Tvars[Tvars %in% vars]
  FN <- Tvars[!(Tvars %in% vars)]
  meanTP <- mean(abs(beta[rownames(beta) %in% TP]))
  meanFN <- mean(abs(beta[rownames(beta) %in% FN]))
  minTP <- min(abs(beta[rownames(beta) %in% TP]))
  minFN <- min(abs(beta[rownames(beta) %in% FN]))
  maxTP <- max(abs(beta[rownames(beta) %in% TP]))
  maxFN <- max(abs(beta[rownames(beta) %in% FN]))
  out <- c(minTP, meanTP, maxTP, minFN, meanFN, maxFN)
  names(out) <- c("minTP", "meanTP", "maxTP", "minFN", "meanFN", "maxFN")
  if (any(is.nan(out))) {
    print("No FN, set NaN and Inf to 0")
    out[is.nan(out)] <- 0
    out[is.infinite(out)] <- 0
  }
  return(out)
}

ReorderCV <- function(varname, vp) {
  order <- which(colSums(vp) != 0)
  order <- order[length(order)]
  for(i in order:1) {
    if (sum(vp[,i] != 0) > 0) {
      row_num <- which(vp[,i] != 0)[1]
      previous <- which(vp[row_num, 1:i-1] != 0)
      if (length(previous) > 0) {
        previous <- previous[length(previous)]
        previous_path <- vp[row_num, previous]
        if (previous_path == 1) {
          order <- c(order, previous)
        } else if (previous_path == 2) {
          order <- c(previous, order)
        }
      } else {
        order <- c(order, i)
      }
    }
  }
  order <- unique(order)
  return(order)
}
# for those that appeared twice but has no kin relationships, maybe do another check on vp
# then, can comment out the order <- c(order,i) and unique（order)
