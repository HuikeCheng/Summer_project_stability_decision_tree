#' Classification And Regression Trees
#'
#' Runs decision trees using implementation from \code{\link[rpart]{rpart}}.
#' This function is not using stability.
#'
#' @inheritParams VariableSelection
#' @param Lambda matrix of parameters controlling the number of splits in the
#'   decision tree.
#' @param ... additional parameters passed to \code{\link[rpart]{rpart}}.
#'
#' @return A list with: \item{selected}{matrix of binary selection status. Rows
#'   correspond to different model parameters. Columns correspond to
#'   predictors.} \item{beta_full}{array of model coefficients. Rows correspond
#'   to different model parameters. Columns correspond to predictors. Indices
#'   along the third dimension correspond to outcome variable(s).}
#'
#' @family underlying algorithm functions
#' @seealso \code{\link{SelectionAlgo}}, \code{\link{VariableSelection}}
#'
#' @references \insertRef{CART}{sharp}
#'
#' @examples
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(pk = 50)
#'
#' # Running the LASSO
#' mycart <- CART(
#'   xdata = simul$xdata,
#'   ydata = simul$ydata,
#'   family = "gaussian"
#' )
#' head(mycart$selected)
#'
#' @export
CART3 <- function(xdata, ydata, Lambda, family, maxsurrogate, ...) {
  # Storing extra arguments
  extra_args <- list(...)
  
  # Defining default Lambda
  Lambda <- cbind(seq(0, 0.2, 0.005))
  
  # Defining the rpart method
  # rpart_method <- switch(family,
  #                        gaussian = "anova",
  #                        binomial = "class",
  #                        cox = "exp",
  #                        poisson = "poisson"
  # )
  
  # Writing the formula
  myformula <- stats::as.formula(paste0("ydata ~ ", paste(paste0("`", colnames(xdata), "`"), collapse = " + ")))
  
  # Creating the data frame
  mydata <- data.frame(ydata = ydata[, 1], xdata)
  
  # Initialising the parameters
  # if (!"cp" %in% names(extra_args)) {
  #   extra_args$cp <- 0
  # }
  # if (!"maxdepth" %in% names(extra_args)) {
  #   extra_args$maxdepth <- 30
  # }
  # if (!any(c("minsplit", "minbucket") %in% names(extra_args))) {
  #   extra_args$minsplit <- 1
  #   extra_args$minbucket <- 1
  # }
  
  # Extracting relevant extra arguments
  # tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = rpart::rpart)
  # tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c(
  #   "formula", "data", "method", "model", "x", "y",
  #   "control", "maxcompete", "maxsurrogate", "xval", "group_x"
  # )]
  
  # Fitting the decision tree
  # mytree <- do.call(partykit::ctree, args = c(
  #   list(
  #     formula = myformula,
  #     data = mydata,
  #    method = rpart_method,
  #    model = FALSE,
  #    x = FALSE,
  #    y = FALSE,
  #    maxcompete = 0,
  #    maxsurrogate = maxsurrogate,
  #    xval = 0
  #   ),
  #   tmp_extra_args
  # ))
  selected <- matrix(0, nrow = nrow(Lambda), ncol = ncol(xdata))
  colnames(selected) <- colnames(xdata)
  
  for (i in 1:nrow(Lambda)) {
    mytree <- do.call(partykit::ctree, args = c(
      list(formula = myformula,
          data = mydata,
          alpha = Lambda[i])))
    
    pvals <- nodeapply(mytree, ids = nodeids(mytree), FUN = function(n) info_node(n)$p.value)
    
    if (list(NULL) %in% pvals) {
      pvals <- Filter(Negate(is.null), pvals)
    }
    
    if (any(pvals < Lambda[i])) {
      selvars <- pvals[which(pvals < Lambda[i])]
      myvars <- sapply(selvars, function(x){names(x[1])})
      selected[i, myvars] <- 1
    }
  }
  
  return(list(selected = selected, beta_full = selected))
}
