SimulateTree <- function(height, n, pk, ev_xy, X=NULL, Y=NULL) {
  ##### Prepare structure of tree
  nlayer <- height - 1
  nrow <- 2^(nlayer)
  ncol <- height
  df <- data.frame(matrix(NA, nrow, ncol))
  colnames(df) <- c(paste0("Level", 1:nlayer), "Leaf")
  
  ##### simulate xdata if not provided
  if(is.null(X)) {
    xdata <- as.data.frame(matrix(rnorm(n*pk),n,pk))
  } else {
    xdata <- as.data.frame(X)
  }
  colnames(xdata) <- paste0("var", 1:pk)
  ##### sample predictors
  # calculate max number of predictors in tree
  nx <- nrow-1
  if(nx > pk) {
    nx <- pk
  }
  # sample nx predictors
  xs <- sample(colnames(xdata), nx)
  # sample predictors in tree from the nx predictors
  for (i in 1:nlayer) {
    df[,i] <- rep(sample(xs, 2^(i-1)), each = 2^(height-i))
  }
  # final set of predictors in tree
  xs <- lapply(df[,-height], unique)
  root <- xs$Level1
  
  ##### simulate subgroup mean
  if (is.null(Y)) {
    df[,height] <- runif(nrow, min = -3, max = 3)
  } else {
    df[,height] <- Y
  }
  
  ##### simulate cutoff value
  uxs <- unique(unlist(xs))
  
  cv <- as.data.frame(matrix(0, length(uxs), nlayer))
  for (i in seq_along(cv)) {
    cv[,i] <- uxs %in% df[,i]
    cv[,i] <- ifelse(cv[,i] == TRUE, 1, NA)
  }
  colnames(cv) <- colnames(df)[-height]
  rownames(cv) <- uxs
  
  for (i in seq_along(uxs)) {
    step <- 1/(sum(!is.na(cv[i,]))+1)
    tmp <- quantile(unlist(xdata[uxs[i]]), probs = seq(0,1,step))
    cv[i,!is.na(cv[i,])] <- sample(tmp[-c(1,length(tmp))], sum(!is.na(cv[i,])))
  }
  
  ##### tree for display
  # df for split
  split <- as.data.frame(matrix(0, nrow, nlayer))
  colnames(split) <- colnames(df)[-height]
  for (i in seq_along(split)) {
    cols <- colnames(split)[i]
    tmp <- df[,i]
    for (j in seq_along(tmp)) {
      varname <- tmp[j]
      split[j,i] <- cv[varname, cols]
    }
  }
  # df for plot
  df_plot <- df
  df_plot$Leaf <- round(df_plot$Leaf, 2)
  for (i in 1:nlayer) {
    df_plot[,i] <- paste(df_plot[,i], ">", round(split[,i], 2), sep = " ")
  }
  
  df_plot$pathString <- rep("", nrow(df_plot))
  for (i in 1:nrow(df_plot)) {
    df_plot$pathString[i] <- paste(unlist(df_plot[i,1:height]), collapse = "/")
  }
  
  tree_display <- as.Node(df_plot)
  
  ##### climbing tree
  # tree for climbing
  df_plot <- df
  df_plot$Leaf <- paste0("Group", 1:(2^nlayer))
  df_plot$pathString <- rep("", nrow(df_plot))
  for (i in 1:nrow(df_plot)) {
    df_plot$pathString[i] <- paste(unlist(df_plot[i,1:height]), collapse = "/")
  }
  
  tree <- as.Node(df_plot)
  
  # set y data for tree
  tree$Set(y = df$Leaf,
           filterFun = isLeaf)
  
  # climb tree to get y data
  tmp <- vector("character", length = n)
  positions <- as.data.frame(matrix(0, n, nlayer))
  colnames(positions) <- paste0("Layer", 1:nlayer)
  
  # layer 1
  col <- tree$Get('name', filterFun = isRoot)
  tmp <- rep(col, n)
  cutoff <- cv[col, 1]
  tmp <- xdata[col] > 0
  tmp <- as.numeric(tmp) + 1 # TRUE -> 2, FALSE -> 1
  positions[,1] <- tmp
  
  # layer 2
  if (nlayer >= 2) {
    for (i in 1:n) {
      col <- tree$Climb(position = tmp[i])$Get('name')[1]
      cutoff <- cv[col, 2]
      tmp[i] <- xdata[i,col] > 0
    }
    tmp <- as.numeric(tmp) + 1
    positions[,2] <- tmp
  }
  
  # layer 3+
  if (nlayer >= 3) {
    for (layer in 2:(nlayer-1)) {
      for (i in 1:n) {
        col <- tree$Climb(position = unlist(positions[i,1:layer]))$Get('name')[1]
        cutoff <- cv[col,layer+1]
        tmp[i] <- xdata[i, col] > 0
      }
      tmp <- as.numeric(tmp) + 1
      positions[,layer+1] <- tmp
    }
  }
  
  # get ypred
  for (i in 1:n) {
    tmp[i] <- tree$Climb(position = unlist(positions[i,1:nlayer]))$y
  }
  
  # generate ydata
  sigma <- sqrt((1 - ev_xy) / ev_xy * stats::var(tmp))
  
  # Sampling from Normal distribution
  ydata <- stats::rnorm(n = n, mean = tmp, sd = sigma)
  ydata <- as.matrix(ydata)
  xdata <- as.matrix(xdata)
  rownames(ydata) <- rownames(xdata) <- paste0("obs", 1:n)
  colnames(ydata) <- "outcome1"
  beta <- rep(0, pk)
  beta[which(beta %in% uxs)] <- 1
  beta <- as.matrix(beta)
  rownames(beta) <- colnames(xdata)
  
  # output
  return(list(xdata = xdata, ydata = ydata, beta = beta, tree = tree_display))
}

# ##### Simulate tree structure
# tree1 <- CreateRegularTree(4, 2)
# SetNodeStyle(tree, style = "rounded", shape = "box")
# plot(tree1)
# print(tree)
# #### Set splits for nodes
# 
# 
# ##### Change to data.frame
# df <- ToDataFrameTypeCol(tree)