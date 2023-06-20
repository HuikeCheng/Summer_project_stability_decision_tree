#' Simulation of undirected graph
#'
#' Simulates the adjacency matrix encoding an unweighted, undirected graph with
#' no self-loops.
#'
#' @inheritParams SimulateGraphical
#' @param pk number of nodes.
#' @param nu expected density (number of edges over the number of node pairs) of
#'   the graph. This argument is only used for \code{topology="random"} or
#'   \code{topology="cluster"} (see argument \code{prob} in
#'   \code{\link[huge]{huge.generator}}).
#' @param ... additional arguments to be passed to
#'   \code{\link[huge]{huge.generator}}.
#'
#' @return A symmetric adjacency matrix encoding an unweighted, undirected graph
#'   with no self-loops.
#'
#' @keywords internal
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


#' Simulation of binary contribution status
#'
#' Simulates the binary contribution status of potential predictor variables
#' from different blocks to outcome variables. For each outcome, the set of true
#' predictors is sampled from one block of potential predictors. If the blocks
#' of variables are independent, the outcomes will be independent too.
#'
#' @inheritParams SimulateSymmetricMatrix
#' @param q number of outcome variables. By default, one block of predictor is
#'   linked to one outcome, i.e. \code{q=sum(pk)}.
#' @param nu vector of probabilities. Each entry corresponds to one block of
#'   predictors and defines the probability for each predictor within the block
#'   to be chosen as true predictor of the corresponding outcome variable.
#' @param orthogonal logical indicating if the outcomes have to be defined from
#'   independent blocks of predictors as encoded in \code{pk}.
#'
#' @return A binary matrix encoding the contribution status of each predictor
#'   variable (columns) to each outcome variable (rows).
#'
#' @keywords internal
SamplePredictors <- function(pk, q = NULL, nu = 0.1, orthogonal = TRUE) {
  # Definition of the number of outcome variables
  if (is.null(q)) {
    q <- length(pk)
  }
  if (length(nu) != q) {
    nu <- rep(nu[1], q)
  }
  
  # Simulation of the binary status for true predictors
  theta <- matrix(0, nrow = q, ncol = sum(pk))
  for (k in 1:q) {
    if (orthogonal) {
      if (k > 1) {
        ids <- seq(cumsum(pk)[k - 1] + 1, cumsum(pk)[k])
      } else {
        ids <- seq(1, cumsum(pk)[k])
      }
      theta[k, ids] <- stats::rbinom(pk[k], size = 1, prob = nu[k])
      
      # Introducing at least one true predictor
      if (sum(theta[k, ids]) == 0) {
        theta[k, sample(ids, size = 1)] <- 1
      }
    } else {
      theta[k, ] <- stats::rbinom(sum(pk), size = 1, prob = nu[k])
      theta[k, k] <- 1
    }
  }
  
  return(t(theta))
}
#' Data simulation for multivariate regression
#'
#' Simulates data with outcome(s) and predictors, where only a subset of the
#' predictors actually contributes to the definition of the outcome(s).
#'
#' @param n number of observations in the simulated dataset. Not used if
#'   \code{xdata} is provided.
#' @param pk number of predictor variables. A subset of these variables
#'   contribute to the outcome definition (see argument \code{nu_xy}). Not used
#'   if \code{xdata} is provided.
#' @param xdata optional data matrix for the predictors with variables as
#'   columns and observations as rows. A subset of these variables contribute to
#'   the outcome definition (see argument \code{nu_xy}).
#' @param family type of regression model. Possible values include
#'   \code{"gaussian"} for continuous outcome(s) or \code{"binomial"} for binary
#'   outcome(s).
#' @param q number of outcome variables.
#' @param theta binary matrix with as many rows as predictors and as many
#'   columns as outcomes. A nonzero entry on row \eqn{i} and column \eqn{j}
#'   indicates that predictor \eqn{i} contributes to the definition of outcome
#'   \eqn{j}.
#' @param nu_xy vector of length \code{q} with expected proportion of predictors
#'   contributing to the definition of each of the \code{q} outcomes.
#' @param beta_abs vector defining the range of nonzero regression coefficients
#'   in absolute values. If \code{continuous=FALSE}, \code{beta_abs} is the set
#'   of possible precision values. If \code{continuous=TRUE}, \code{beta_abs} is
#'   the range of possible precision values. Note that regression coefficients
#'   are re-scaled if \code{family="binomial"} to ensure that the desired
#'   concordance statistic can be achieved (see argument \code{ev_xy}) so they
#'   may not be in this range.
#' @param beta_sign vector of possible signs for regression coefficients.
#'   Possible inputs are: \code{1} for positive coefficients, \code{-1} for
#'   negative coefficients, or \code{c(-1, 1)} for both positive and negative
#'   coefficients.
#' @param continuous logical indicating whether to sample regression
#'   coefficients from a uniform distribution between the minimum and maximum
#'   values in \code{beta_abs} (if \code{continuous=TRUE}) or from proposed
#'   values in \code{beta_abs} (if \code{continuous=FALSE}).
#' @param ev_xy vector of length \code{q} with expected goodness of fit measures
#'   for each of the \code{q} outcomes. If \code{family="gaussian"}, the vector
#'   contains expected proportions of variance in each of the \code{q} outcomes
#'   that can be explained by the predictors. If \code{family="binomial"}, the
#'   vector contains expected concordance statistics (i.e. area under the ROC
#'   curve) given the true probabilities.
#'
#' @return A list with: \item{xdata}{input or simulated predictor data.}
#'   \item{ydata}{simulated outcome data.} \item{beta}{matrix of true beta
#'   coefficients used to generate outcomes in \code{ydata} from predictors in
#'   \code{xdata}.} \item{theta}{binary matrix indicating the predictors from
#'   \code{xdata} contributing to the definition of each of the outcome
#'   variables in \code{ydata}.}
#'
#' @family simulation functions
#'
#' @references \insertRef{ourstabilityselection}{fake}
#'
#' @examples
#' \donttest{
#' ## Independent predictors
#'
#' # Univariate continuous outcome
#' set.seed(1)
#' simul <- SimulateRegression(pk = 15)
#' summary(simul)
#'
#' # Univariate binary outcome
#' set.seed(1)
#' simul <- SimulateRegression(pk = 15, family = "binomial")
#' table(simul$ydata)
#'
#' # Multiple continuous outcomes
#' set.seed(1)
#' simul <- SimulateRegression(pk = 15, q = 3)
#' summary(simul)
#'
#'
#' ## Blocks of correlated predictors
#'
#' # Simulation of predictor data
#' set.seed(1)
#' xsimul <- SimulateGraphical(pk = rep(5, 3), nu_within = 0.8, nu_between = 0, v_sign = -1)
#' Heatmap(cor(xsimul$data),
#'   legend_range = c(-1, 1),
#'   col = c("navy", "white", "darkred")
#' )
#'
#' # Simulation of outcome data
#' simul <- SimulateRegression(xdata = xsimul$data)
#' print(simul)
#' summary(simul)
#'
#'
#' ## Choosing expected proportion of explained variance
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 1000, pk = 15, q = 3, ev_xy = c(0.9, 0.5, 0.2))
#' summary(simul)
#'
#' # Comparing with estimated proportion of explained variance
#' summary(lm(simul$ydata[, 1] ~ simul$xdata))
#' summary(lm(simul$ydata[, 2] ~ simul$xdata))
#' summary(lm(simul$ydata[, 3] ~ simul$xdata))
#'
#'
#' ## Choosing expected concordance (AUC)
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(
#'   n = 500, pk = 10,
#'   family = "binomial", ev_xy = 0.9
#' )
#'
#' # Comparing with estimated concordance
#' fitted <- glm(simul$ydata ~ simul$xdata,
#'   family = "binomial"
#' )$fitted.values
#' Concordance(observed = simul$ydata, predicted = fitted)
#' }
#' @export
SimulateInter <- function(n = 100, pk = 10, xdata = NULL,
                              family = "gaussian", q = 1,
                              theta = NULL, nu_xy = 0.2,
                              beta_abs = c(0.1, 1), beta_sign = c(-1, 1), continuous = TRUE,
                              ev_xy = 0.7,
                              type = "Pair",
                              Number = 1,
                              proportion = 1) {
  # TODO in future versions: introduce more families ("multinomial" and "cox")
  # Checking that either n and pk or xdata are provided
  if (is.null(xdata)) {
    if (is.null(pk) | is.null(n)) {
      stop("Argument 'xdata' must be provided if 'pk' and 'n' are not provided.")
    }
  }
  
  # Checking other inputs
  if (length(ev_xy) != q) {
    ev_xy <- rep(ev_xy[1], q)
  }
  if (length(nu_xy) != q) {
    nu_xy <- rep(nu_xy[1], q)
  }
  
  # Creating objects not provided as input
  if (!is.null(xdata)) {
    n <- nrow(xdata)
    p <- ncol(xdata)
  } else {
    p <- sum(pk)
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
  }
  
  # Checking theta if provided
  if (!is.null(theta)) {
    if (q == 1) {
      if (is.vector(theta)) {
        theta <- cbind(theta)
      }
    }
    if (ncol(theta) != q) {
      stop("Arguments 'theta' and 'q' are not compatible. Please provide a matrix 'theta' with 'q' columns.")
    }
    if (nrow(theta) != p) {
      stop("Please provide a matrix 'theta' with as many columns as predictors.")
    }
    theta <- ifelse(theta != 0, yes = 1, no = 0)
  }
  
  # Sampling true predictors
  if (is.null(theta)) {
    theta <- SamplePredictors(pk = p, q = q, nu = nu_xy, orthogonal = FALSE)
  }
  
  # Sampling regression coefficients
  beta <- theta
  if (continuous) {
    beta <- beta * matrix(stats::runif(n = nrow(beta) * ncol(beta), min = min(beta_abs), max = max(beta_abs)),
                          nrow = nrow(beta), ncol = ncol(beta)
    )
  } else {
    beta <- beta * matrix(base::sample(beta_abs, size = nrow(beta) * ncol(beta), replace = TRUE),
                          nrow = nrow(beta), ncol = ncol(beta)
    )
  }
  beta <- beta * matrix(base::sample(beta_sign, size = nrow(beta) * ncol(beta), replace = TRUE),
                        nrow = nrow(beta), ncol = ncol(beta)
  )
  
  # Changing the relationship between true predictor and y
  if (q == 1) {
    id <- which(theta == 1)
    newx <- vector("list", length = Number)
    j <- 1
    if (type == "Pair") {
      for (i in 1:Number) {
        newx[[i]] <- xdata[,id[j]] * xdata[,id[j+1]]
        j+2
      }
      newx <- list.cbind(newx)
      xdata <- cbind(xdata,newx)
      colnames(xdata) <- paste0("var", 1:(pk+Number))
      theta <- as.matrix(c(theta, rep(1,Number)))
      beta <- as.matrix(c(beta, rep(1,Number)))
    } else if (type == "Layer") {
      id <- sample(id, Number)
      newx <- rep(1, n)
      for (i in 1:Number) {
        newx <- newx * xdata[,id[i]]
      }
      newx <- as.matrix(newx)
      xdata <- cbind(xdata, newx)
      colnames(xdata) <- paste0("var", 1:(pk+1))
      theta <- as.matrix(c(theta, 1))
      beta <- as.matrix(c(beta, 1))
    }
  } else {
    stop("More than 1 outcome")
  }
  # Sampling outcome data
  ydata <- matrix(NA, ncol = q, nrow = nrow(xdata))
  if (family == "gaussian") {
    for (j in 1:q) {
      # Linear combination
      ypred <- xdata %*% beta[, j]
      
      # Calculating standard deviation that achieves desired proportion of explained variance
      sigma <- sqrt((1 - ev_xy[j]) / ev_xy[j] * stats::var(ypred)) # using var(ypred) is a shortcut (could use true variances but not always available)
      
      # Sampling from Normal distribution
      ydata[, j] <- stats::rnorm(n = n, mean = ypred, sd = sigma)
    }
  }
  if (family == "binomial") {
    for (j in 1:q) {
      # Linear combination
      crude_log_odds <- xdata %*% beta[, j]
      
      # Identifying a relevant range of scaling factors (logistic distribution)
      s_max <- max(abs(crude_log_odds)) / log(0.51 / 0.49) # scale that gives all probabilities between 0.49 and 0.51 (expected c stat close to 0.5)
      s_min <- min(abs(crude_log_odds)) / log(0.99 / 0.01) # scale that gives all probabilities above 0.99 or below 0.01 (expected c stat close to 1)
      
      # Finding scaling factor that gives desired AUC (interval required to ensure optimisation works)
      argmax_scaling_factor <- stats::optimise(
        f = TuneCStatisticLogit,
        crude_log_odds = crude_log_odds,
        auc = ev_xy[j],
        lower = 1 / s_max, upper = 1 / s_min
      )
      scaling_factor <- argmax_scaling_factor$minimum
      
      # Applying scaling factor
      beta[, j] <- beta[, j] * scaling_factor
      log_odds <- crude_log_odds * scaling_factor
      
      # Calculating probabilities from log-odds (inverse logit)
      proba <- 1 / (1 + exp(-log_odds))
      
      # Sampling from Bernouilli distribution
      ydata[, j] <- stats::rbinom(n = n, size = 1, prob = proba)
    }
  }
  
  # Defining row and column names of ydata
  rownames(ydata) <- rownames(xdata)
  colnames(ydata) <- paste0("outcome", 1:q)
  
  # Defining row and column names of beta and theta
  rownames(beta) <- rownames(theta) <- colnames(xdata)
  colnames(beta) <- colnames(theta) <- colnames(ydata)
  
  # Preparing the output
  out <- list(xdata = xdata, ydata = ydata, beta = beta, theta = theta)
  
  # Defining the class
  class(out) <- "simulation_regression"
  
  return(out)
}
