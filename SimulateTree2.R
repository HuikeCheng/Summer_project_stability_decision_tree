SimulateTree <- function(height, n, pk, ev_xy, X=NULL, Y=NULL, Y_abs = 3, repeated_pred = TRUE, multivariate_normal = TRUE) {
  ##### Prepare structure of tree as dataframe
  nlayer <- height - 1
  ngroups <- 2^(nlayer)
  df <- data.frame(matrix(NA, ngroups, height))
  colnames(df) <- c(paste0("Level", 1:nlayer), "Leaf")
  
  ##### simulate xdata if not provided
  if(!is.null(X)) {
    xdata <- as.matrix(X)
  } else {
    if (multivariate_normal == TRUE) {
      xsimul <- SimulateGraphical(
        n = n, pk = pk, theta = NULL,
        implementation = HugeAdjacency, topology = "random",
        nu_within = 0, nu_between = 0, nu_mat = NULL,
        v_within = 0, v_between = 0,
        v_sign = c(-1, 1), continuous = TRUE,
        pd_strategy = "diagonally_dominant", ev_xx = NULL, scale_ev = TRUE,
        u_list = c(1e-10, 1), tol = .Machine$double.eps^0.25,
        scale = TRUE, output_matrices = FALSE
      )
      xdata <- xsimul$data
      print("mn")
    } else {
      xdata <- matrix(rnorm(n*pk),n,pk)
      colnames(xdata) <- paste0("var", 1:pk)
    }
  }
  
  ##### sample predictors
  if (repeated_pred == TRUE) {
    # calculate max number of predictors in tree
    mn <- ngroups-1
    if(mn > ncol(xdata)) {
      mn <- ncol(xdata)
    }
    # sample mn number of candidate predictors
    xs <- sample(1:ncol(xdata), mn)
    # sample predictors in tree from the mn predictors
    for (i in 1:nlayer) {
      df[,i] <- rep(sample(xs, 2^(i-1)), each = 2^(height-i))
    }
  } else {
    # sample predictors in tree from colnames(xdata)
    for (i in 1:nlayer) {
      df[,i] <- rep(sample(1:ncol(xdata), 2^(i-1)), each = 2^(height-i))
    }
  }
  # predictors in tree
  xs <- lapply(df[,-height], unique)
  
  ##### simulate subgroup mean
  if (!is.null(Y)) {
    if (length(Y) == ngroups) {
      df[,height] <- Y
    } else {
      stop("number of subgroup means provided does not equal number of subgroups")
    }
  } else {
    df[,height] <- runif(ngroups, min = -Y_abs, max = Y_abs)
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
    tmp <- quantile(unlist(xdata[,uxs[i]]), probs = seq(0,1,step))
    cv[i,!is.na(cv[i,])] <- sample(tmp[-c(1,length(tmp))], sum(!is.na(cv[i,])))
  }
  
  # reorder cutoff values
  structure <- df[!duplicated(df[,1:nlayer]),1:nlayer]
  path <- structure
  for (i in 1:(nlayer-1)) {
    path[,i] <- rep(rep(c(1,2), each = 2^(nlayer-1-i)), 2^(i-1))
  }
  path[,nlayer] <- rep(3, nrow(path))
  
  for (i in seq_along(uxs)) {
    if (sum(!is.na(cv[i,])) > 1) {
      vp <- apply(structure, 2, FUN = function(x) {ifelse(x != uxs[i], 0, 1)})
      vp <- vp*path
      if (any(apply(vp, 1, FUN = function(x){sum(x!=0)}) > 1)) {
        order <- reorderCV(vp = vp)
        values <- unlist(cv[i, !is.na(cv[i,])])
        cv[i, order] <- sort(values)
      }
    }
  }
  
  ##### tree for display
  # df for split
  split <- as.data.frame(matrix(0, ngroups, nlayer))
  colnames(split) <- colnames(df)[-height]
  for (i in seq_along(split)) {
    tmp <- as.character(df[,i])
    split[,i] <- cv[tmp,i]
  }
  # df for plot
  df_plot <- df
  df_plot$Leaf <- round(df_plot$Leaf, 2)
  for (i in 1:nlayer) {
    df_plot[,i] <- paste("var", df_plot[,i], " > ", round(split[,i], 2), sep = "")
  }
  
  ##### climbing tree
  # climb tree to get y data
  tmp <- vector("character", length = n)
  positions <- as.data.frame(matrix(NA, n, nlayer))
  colnames(positions) <- paste0("Level", 1:nlayer)
  
  # layer 1
  col <- xs[[1]] # var at root
  cutoff <- cv[as.character(col), 1]
  tmp <- xdata[,col] > cutoff
  positions[,1] <- as.numeric(tmp) # TRUE -> 1, FALSE -> 0
  
  # layer 2
  index <- positions[,1] + 1
  col <- xs[[2]][index]
  cutoff <- cv[as.character(col),2]
  tmp <- xdata[cbind(1:n,col)] > cutoff
  positions[,2] <- as.numeric(tmp)
  
  # layer 3+
  for (i in 3:nlayer) {
    index <- apply(positions[,1:(i-1)], 1, Bin2Dec) + 1
    col <- xs[[i]][index]
    cutoff <- cv[as.character(col),i]
    tmp <- xdata[cbind(1:n, col)] > cutoff
    positions[,i] <- as.numeric(tmp)
  }
  
  # get ypred
  groups <- apply(positions[,1:nlayer], 1, Bin2Dec) + 1
  ypred <- df[groups,height]
  
  # generate ydata
  sigma <- sqrt((1 - ev_xy) / ev_xy * stats::var(ypred))
  
  # Sampling from Normal distribution
  ydata <- stats::rnorm(n = n, mean = ypred, sd = sigma)
  ydata <- as.matrix(ydata)
  rownames(ydata) <- rownames(xdata) <- paste0("obs", 1:n)
  colnames(ydata) <- "outcome1"
  
  ##### check number of obs in each subgroup
  size <- summary(factor(groups))
  if (sum(!(1:ngroups) %in% names(size)) > 0) {
    # empty groups
    empty_groups <- which(!(1:ngroups) %in% names(size))
    df$emptycol <- rep("", ngroups)
    df_plot$emptycol <- rep("", ngroups)
    
    # pruning
    df <- df[-empty_groups,]
    df_plot <- df_plot[-empty_groups,]
    
    # remove extraneous split
    for (layer in nlayer:1) {
      for (i in unique(df[,layer])) {
        index <- which(df[,layer] == i)
        numchild <- length(unique(df[index,(layer+1)]))
        if (numchild < 2) {
          df[index,layer:height] <- df[index,(layer+1):(height+1)]
          df_plot[index,layer:height] <- df_plot[index,(layer+1):(height+1)]
        }
      }
    }
    # remove emptycol used for pruning
    df <- df[,-(height+1)]
    df_plot <- df_plot[,-(height+1)]
    # check if final height is still height asked
    emptylevels <- which(apply(df,2,FUN = function(x){sum(x!="")}) == 0)
    if (length(emptylevels) > 0) {
      df <- df[,-emptylevels]
      df_plot <- df_plot[,-emptylevels]
      warning(paste0("Final tree height = ", height-(length(emptylevels))))
    }
  }
  
  # get tree for display
  df_plot$pathString <- apply(df_plot, 1, FUN=function(x){paste0(x,collapse="/")})
  tree_display <- as.Node(df_plot)
  SetNodeStyle(tree_display, style = "rounded", shape = "box")
  
  ##### set theta
  df$pathString <- apply(df, 1, FUN=function(x){paste0(x,collapse="/")})
  tree <- as.Node(df)
  var_final <- as.numeric(unique(tree$Get('name', filterFun = isNotLeaf)))
  theta <- rep(0, pk)
  names(theta) <- colnames(xdata)
  theta[var_final] <- 1
  
  # output
  return(list(xdata = xdata, ydata = ydata, theta = theta, tree = tree_display))
}

#####################
# library(data.tree)
# library(fake)
# source("Functions.R")
# 
# height <- 5
# n <- 1000
# pk <- 500
# ev_xy <- 0.1
# 
# start.time <- Sys.time()
# set.seed(4)
# simul <- SimulateTree(height, n, pk, 0.2)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# # 4:  s
# # 5: 7.42 s
# plot(simul$tree)

# things to work on
# add parameter to adjust allow repeat or not. done
# check subgroup means, if too close, not meaningful? weak predictor
# another check, when sample size too small, height to big, all groups could be empty? added warning for height
# check at the end if tree has any problems, if data correctly end up in groups, and if theta corresponds to tree

# 
# ##### Change to data.frame
# df <- ToDataFrameTypeCol(tree)