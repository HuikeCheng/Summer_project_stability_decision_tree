SimulateTree <- function(height, n, pk, ev_xy, X=NULL, Y=NULL) {
  ##### Prepare structure of tree as dataframe
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
  mn <- nrow-1
  if(mn > pk) {
    mn <- pk
  }
  # sample mn number of candidate predictors
  xs <- sample(colnames(xdata), mn) # add argument to allow repetitions or not
  # sample predictors in tree from the mn predictors
  for (i in 1:nlayer) {
    df[,i] <- rep(sample(xs, 2^(i-1)), each = 2^(height-i))
  }
  # final set of predictors in tree
  xs <- lapply(df[,-height], unique)
  
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
  
  # sample cutoff value from quantiles
  for (i in seq_along(uxs)) {
    step <- 1/(sum(!is.na(cv[i,]))+1)
    tmp <- quantile(unlist(xdata[uxs[i]]), probs = seq(0,1,step))
    cv[i,!is.na(cv[i,])] <- sample(tmp[-c(1,length(tmp))], sum(!is.na(cv[i,])))
  }
  
  # reorder cutoff values
  structure <- df[!duplicated(df[,1:nlayer]),1:nlayer]
  path <- structure
  for (i in 1:(nlayer-1)) {
    path[,i] <- rep(rep(c(1,2), each = 2^(nlayer-1-i)), 2^(i-1))
  }
  path[,nlayer] <- rep(3, nrow(path))
  
  for (i in uxs) {
    if (sum(!is.na(cv[i,])) > 1) {
      vp <- apply(structure, 2, FUN = function(x) {ifelse(x != i, 0, 1)})
      vp <- vp*path
      if (any(apply(vp, 1, FUN = function(x){sum(x!=0)})) > 1) {
        order <- reorderCV(varname = i, vp = vp)
        values <- unlist(cv[i, !is.na(cv[i,])])
        cv[i, order] <- sort(values)
      }
    }
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
  df_tree <- df
  df_tree$Leaf <- paste0("Group", 1:(2^nlayer))
  df_tree$pathString <- rep("", nrow(df_tree))
  for (i in 1:nrow(df_tree)) {
    df_tree$pathString[i] <- paste(unlist(df_tree[i,1:height]), collapse = "/")
  }
  
  tree <- as.Node(df_tree)
  
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
  tmp <- xdata[col] > cutoff
  tmp <- as.numeric(tmp) + 1 # TRUE -> 2, FALSE -> 1
  positions[,1] <- tmp
  
  # layer 2
  if (nlayer >= 2) {
    for (i in 1:n) {
      col <- tree$Climb(position = tmp[i])$Get('name')[1]
      cutoff <- cv[col, 2]
      tmp[i] <- as.numeric(xdata[i,col] > cutoff)
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
        tmp[i] <- as.numeric(xdata[i, col] > cutoff)
      }
      tmp <- as.numeric(tmp) + 1
      positions[,layer+1] <- tmp
    }
  }
  
  # get ypred
  groups <- vector("numeric", length = n)
  for (i in 1:n) {
    tmp[i] <- tree$Climb(position = unlist(positions[i,1:nlayer]))$y
    groups[i] <- tree$Climb(position = unlist(positions[i,1:nlayer]))$Get('name')
  }
  
  # generate ydata
  sigma <- sqrt((1 - ev_xy) / ev_xy * stats::var(tmp))
  
  # Sampling from Normal distribution
  ydata <- stats::rnorm(n = n, mean = tmp, sd = sigma)
  ydata <- as.matrix(ydata)
  xdata <- as.matrix(xdata)
  rownames(ydata) <- rownames(xdata) <- paste0("obs", 1:n)
  colnames(ydata) <- "outcome1"
  
  ##### check number of obs in each subgroup
  size <- summary(factor(groups))
  if (sum(!df_tree$Leaf %in% names(size)) > 0) {
    empty_group <- df_tree$Leaf[which(!df_tree$Leaf %in% names(size))]
    empty_leaf <- df_plot$Leaf[which(!df_tree$Leaf %in% names(size))]
    for (i in seq_along(empty_group)) {
      Prune(tree, function(x) {!(empty_group[i] %in% x$Get('name') & x$level == nlayer)})
      Prune(tree_display, function(x) {!(empty_leaf[i] %in% x$Get('name') & x$level == nlayer)})
    }
  }
  var_final <- unique(tree$Get('name', filterFun = isNotLeaf))
  
  ##### set theta
  theta <- rep(0, pk)
  theta <- as.matrix(theta)
  rownames(theta) <- colnames(xdata)
  theta[which(rownames(theta) %in% var_final)] <- 1
  
  # output
  return(list(xdata = xdata, ydata = ydata, theta = theta, tree = tree_display))
}

#
height <- 5
n <- 1000
pk <- 500
ev_xy <- 0.1
X <- NULL
Y <- NULL

start.time <- Sys.time()
set.seed(1)
simul <- SimulateTree(height, n, pk, ev_xy, X, Y)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken 
# 4: 6.12 s
# 5: 5.95 s
plot(simul$tree)

# things to work on further adjustments, or alternatives to pruning
# add parameter to adjust allow repeat or not
# check subgroup means, if too close, not meaningful?
# another check, when sample size too small, height to big, all groups could be empty?


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

library(data.tree)

df <- data.frame(matrix(NA, nrow=4,ncol=3))
df[,1] <- rep("node1", 4)
df[,2] <- rep(c("node2", "node3"), each=2)
df[,3] <- c("node4", "node5", "node6", "node7")

df$pathString <- apply(df, 1, FUN=function(x){paste0(x,collapse="/")})

tree <- as.Node(df)
SetNodeStyle(tree, style = "rounded", shape = "box")
plot(tree)

Prune(tree, pruneFun = function(x) {!("node4" %in% x$Get("name") & isLeaf(x))})
plot(tree)

Prune(tree, pruneFun = function(x) {!("node4" %in% x$Get("name") & x$level == 2)})
plot(tree)