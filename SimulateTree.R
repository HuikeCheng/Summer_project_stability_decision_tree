SimulateTree <- function(height, n, pk, ev_xy) {
  ##### Build a tree
  nlayer <- height - 1
  nrow <- 2^(nlayer)
  ncol <- height
  df <- data.frame(matrix(NA, nrow, ncol))
  colnames(df) <- c(paste0("Level", 1:nlayer), "Leaf")
  df[,1] <- rep("var1", nrow)
  df[,2] <- rep(c("var2", "var3"), each = 4)
  df[,3] <- rep(c("var4", "var5", "var6", "var7"), each = 2)
  df[,4] <- runif(nrow, min = -3, max = 3)
  
  ##### split
  split <- data.frame(matrix(0, nrow, nlayer))
  
  ####### tree for show
  df_plot <- df
  df_plot[,1] <- paste(df_plot[,1], ">", split[,1], sep = " ")
  df_plot[,2] <- paste(df_plot[,2], ">", split[,2], sep = " ")
  df_plot[,3] <- paste(df_plot[,3], ">", split[,3], sep = " ")
  df_plot$Leaf <- round(df_plot$Leaf, 2)
  df_plot$pathString <- paste(df_plot[1,1], 
                              df_plot[,2], 
                              df_plot[,3],
                              df_plot[,4],
                              sep = "/")
  tree_display <- as.Node(df_plot)
  
  ########## tree for climbing
  df_plot <- df
  df_plot$Leaf <- round(df_plot$Leaf, 2)
  df_plot$pathString <- paste(df_plot[,1], 
                              df_plot[,2], 
                              df_plot[,3],
                              df_plot[,4],
                              sep = "/")
  tree <- as.Node(df_plot)
  # set y data for tree
  tree$Set(y = df$Leaf,
           filterFun = isLeaf)
  
  # xdata
  X <- as.data.frame(matrix(rnorm(n*pk),n,pk))
  colnames(X) <- paste0("var", 1:ncol(X))
  
  # climb tree to get y data
  tmp <- vector("character", length = n)
  positions <- as.data.frame(matrix(0, n, height-1))
  colnames(positions) <- paste0("Layer", 1:(height-1))
  
  # layer 1
  col <- tree$Get('name', filterFun = isRoot)
  tmp <- rep(col, 1000)
  tmp <- X[col] <= 0
  tmp <- as.numeric(tmp) + 1 # TRUE -> 2, FALSE -> 1
  positions[,1] <- tmp
  
  # layer 2
  for (i in 1:n) {
    col <- tree$Climb(position = tmp[i])$Get('name')[1]
    tmp[i] <- X[i,col] <= 0
  }
  tmp <- as.numeric(tmp) + 1
  positions[,2] <- tmp
  
  # layer 3
  for (i in 1:n) {
    col <- tree$Climb(position = unlist(positions[i,1:2]))$Get('name')[1]
    tmp[i] <- X[i, col] <= 0
  }
  tmp <- as.numeric(tmp) + 1
  positions[,3] <- tmp
  
  # get ypred
  for (i in 1:n) {
    tmp[i] <- tree$Climb(position = unlist(positions[i,1:3]))$y
  }
  
  # generate ydata
  sigma <- sqrt((1 - ev_xy) / ev_xy * stats::var(tmp)) # using var(ypred) is a shortcut (could use true variances but not always available)
  
  # Sampling from Normal distribution
  ydata <- stats::rnorm(n = n, mean = tmp, sd = sigma)
  ydata <- as.matrix(ydata)
  xdata <- as.matrix(X)
  rownames(ydata) <- rownames(xdata) <- paste0("obs", 1:n)
  colnames(ydata) <- "outcome1"
  beta <- rep(0, pk)
  beta[1:7] <- 1
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