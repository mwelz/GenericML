#' Calculate lower and upper median
#'
#' Calculates the lower and and median of a vector as proposed in Comment 4.2 in the paper.
#'
#' @param x A numeric vector.
#'
#' @return
#' A list with the upper and lower median and the Med statistic (which is their mean).
#'
#' @references
#' Chernozhukov V., Demirer M., Duflo E., Fernández-Val I. (2020). \dQuote{Generic Machine Learning Inference on Heterogenous Treatment Effects in Randomized Experiments.} \emph{arXiv preprint arXiv:1712.04802}. URL: \url{https://arxiv.org/abs/1712.04802}.
#'
#' @examples
#' set.seed(1)
#' x <- runif(100)
#' Med(x)
#'
#' @export
Med <- function(x){

  # get the empirical CDF of X
  ecdf.x <- stats::ecdf(x)

  # evaluate the ensuing probabilities
  F.x <- ecdf.x(x)

  # get lower median and upper median
  lower       <- min(x[F.x >= 0.5])
  upper.array <- x[(1 - F.x) >= 0.5]
  upper       <- ifelse(length(upper.array) == 0,
                        lower,
                        max(upper.array)) # account for case where upper.array is empty

  return(list(lower_median = lower,
              upper_median = upper,
              Med = mean(c(lower, upper))))

} # END FUN



#' Partition a vector into quantile groups
#'
#' Partitions a vector into quantile groups and returns a logical matrix indicating group membership.
#'
#' @param x A numeric vector to be partitioned.
#' @param cutoffs A numeric vector denoting the quantile cutoffs for the partition. Default are the quartiles: \code{c(0.25, 0.5, 0.75)}.
#'
#' @return
#' An object of type \code{"quantile_group"}, which is a logical matrix indicating group membership.
#'
#' @examples
#' set.seed(1)
#' x <- runif(100)
#' cutoffs <- c(0.25, 0.5, 0.75)
#' quantile_group(x, cutoffs)
#'
#' @export
quantile_group <- function(x,
                           cutoffs = c(0.25, 0.5, 0.75)){


  # input checks
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(cutoffs))
  stopifnot(0 < min(cutoffs) & max(cutoffs) < 1)

  # return
  quantile_group_NoChecks(x = x, cutoffs = cutoffs)


} # FUN


# same as above, just w/o input checks
quantile_group_NoChecks <- function(x = x,
                                    cutoffs = cutoffs){

  # get quantiles
  q         <- stats::quantile(x, cutoffs)
  q         <- c(-Inf, q, Inf)

  # check if breaks are unique: if x exhibits low variation, there might be empty quantile bins, which can cause an error in the cut() function. In this case, we add random noise to x to induce variation. NB: this bug has been spotted and fixed by Lucas Kitzmueller. All credits for this fix go to him!
  if(length(unique(q)) != length(q)){
    # specify standard deviation of the noise (x may have zero variation)
    sd <- ifelse(stats::var(x) == 0, 0.001, sqrt(stats::var(x) / 20))
    # add noise and updare quantiles
    x <- x + stats::rnorm(length(x), mean = 0, sd = sd)
    q <- stats::quantile(x, cutoffs)
    q <- c(-Inf, q, Inf)
  } # IF

  groups    <- as.character(cut(x, breaks = q, include.lowest = TRUE, right = FALSE, dig.lab = 3))
  group.nam <- unique(groups)
  group.nam <- group.nam[order(
    as.numeric(substr(sub("\\,.*", "", group.nam), 2, stop = 1e8L)),
    decreasing = FALSE)] # ensure the order is correct

  # get the grouping matrix
  group.mat <- sapply(1:length(group.nam), function(j) groups == group.nam[j])
  colnames(group.mat) <- gsub(",", ", ", gsub(" ", "", group.nam))

  # return
  return(structure(group.mat, type = "quantile_group"))

} # FUN


#' Single iteration of the GenericML algorithm
#'
#' Performs generic ML inference for a single learning technique and a given split of the data. Can be seen as a single iteration of Algorithm 1 in the paper.
#'
#' @param Z A numeric design matrix that holds the covariates in its columns.
#' @param D A binary vector of treatment assignment. Value one denotes assignment to the treatment group and value zero assignment to the control group.
#' @param Y A numeric vector containing the response variable.
#' @param learner A character specifying the machine learner to be used for estimating the baseline conditional average (BCA) and conditional average treatment effect (CATE). Either \code{'lasso'}, \code{'random_forest'}, \code{'tree'}, or a custom learner specified with \code{mlr3} syntax. In the latter case, do \emph{not} specify in the \code{mlr3} syntax specification if the learner is a regression learner or classification learner. Example: \code{'mlr3::lrn("ranger", num.trees = 100)'} for a random forest learner with 100 trees. Note that this is a string and the absence of the \code{classif.} or \code{regr.} keywords. See \url{https://mlr3learners.mlr-org.com} for a list of \code{mlr3} learners.
#' @param propensity_scores A numeric vector of propensity score estimates.
#' @param M_set A numerical vector of indices of observations in the main sample.
#' @param A_set A numerical vector of indices of observations in the auxiliary sample. Default is complementary set to \code{M_set}.
#' @param Z_CLAN A numeric matrix holding variables on which classification analysis (CLAN) shall be performed. CLAN will be performed on each column of the matrix. If \code{NULL} (default), then \code{Z_CLAN = Z}, i.e. CLAN is performed for all variables in \code{Z}.
#' @param HT Logical. If \code{TRUE}, a Horvitz-Thompson (HT) transformation is applied in the BLP and GATES regressions. Default is \code{FALSE}.
#' @param quantile_cutoffs The cutoff points of the quantiles that shall be used for GATES grouping. Default is \code{c(0.25, 0.5, 0.75)}, which corresponds to the four quartiles.
#' @param X1_BLP Specifies the design matrix \eqn{X_1} in the regression. Must be an object of class  \code{"\link{setup_X1}"}. See the documentation of \code{\link{setup_X1}()} for details.
#' @param X1_GATES Same as \code{X1_BLP}, just for the GATES regression.
#' @param diff_GATES Specifies the generic targets of GATES. Must be an object of class \code{"\link{setup_diff}"}. See the documentation of \code{\link{setup_diff}()} for details.
#' @param diff_CLAN Same as \code{diff_GATES}, just for the CLAN generic targets.
#' @param vcov_BLP Specifies the covariance matrix estimator in the BLP regression. Must be an object of class \code{"\link{setup_vcov}"}. See the documentation of \code{\link{setup_vcov}()} for details.
#' @param vcov_GATES Same as \code{vcov_BLP}, just for the GATES regression.
#' @param equal_variances_CLAN Logical. If \code{TRUE}, then all within-group variances of the CLAN groups are assumed to be equal. Default is \code{FALSE}. This specification is required for heteroskedasticity-robust variance estimation on the difference of two CLAN generic targets (i.e. variance of the difference of two means). If \code{TRUE} (corresponds to homoskedasticity assumption), the pooled variance is used. If \code{FALSE} (heteroskedasticity), the variance of Welch's t-test is used.
#' @param significance_level Significance level for VEIN. Default is 0.05.
#' @param min_variation Specifies a threshold for the minimum variation of the BCA/CATE predictions. If the variation of a BCA/CATE prediction falls below this threshold, random noise with distribution \eqn{N(0, var(Y)/20)} is added to it. Default is \code{1e-05}.
#'
#' @return
#' A list with the following components:
#' \describe{
#' \item{\code{BLP}}{An object of class \code{"\link{BLP}"}.}
#' \item{\code{GATES}}{An object of class \code{"\link{GATES}"}.}
#' \item{\code{CLAN}}{An object of class \code{"\link{CLAN}"}.}
#' \item{\code{proxy_BCA}}{An object of class \code{"\link{proxy_BCA}"}.}
#' \item{\code{proxy_CATE}}{An object of class \code{"\link{proxy_CATE}"}.}
#' \item{\code{best}}{Estimates of the \eqn{\Lambda} parameters for finding the best learner. Returned by \code{\link{lambda_parameters}()}.}
#' }
#'
#' @details
#' The specifications \code{"lasso"}, \code{"random_forest"}, and \code{"tree"} in \code{learner} correspond to the following \code{mlr3} specifications (we omit the keywords \code{classif.} and \code{regr.}). \code{"lasso"} is a cross-validated Lasso estimator, which corresponds to \code{'mlr3::lrn("cv_glmnet", s = "lambda.min", alpha = 1)'}. \code{"random_forest"} is a random forest with 500 trees, which corresponds to \code{'mlr3::lrn("ranger", num.trees = 500)'}. \code{"tree"} is a tree learner, which corresponds to \code{'mlr3::lrn("rpart")'}.
#'
#' @references
#' Chernozhukov V., Demirer M., Duflo E., Fernández-Val I. (2020). \dQuote{Generic Machine Learning Inference on Heterogenous Treatment Effects in Randomized Experiments.} \emph{arXiv preprint arXiv:1712.04802}. URL: \url{https://arxiv.org/abs/1712.04802}.
#'
#' Lang M., Binder M., Richter J., Schratz P., Pfisterer F., Coors S., Au Q., Casalicchio G., Kotthoff L., Bischl B. (2019). \dQuote{mlr3: A Modern Object-Oriented Machine Learning Framework in R.} \emph{Journal of Open Source Software}, \bold{4}(44), 1903. \doi{10.21105/joss.01903}.
#'
#' @seealso
#' \code{\link{GenericML}()}
#'
#' @examples
#' if(require("ranger")){
#' ## generate data
#' set.seed(1)
#' n  <- 150                        # number of observations
#' p  <- 5                          # number of covariates
#' Z  <- matrix(runif(n*p), n, p)   # design matrix
#' D  <- rbinom(n, 1, 0.5)          # random treatment assignment
#' Y  <- runif(n)                   # outcome variable
#' propensity_scores <- rep(0.5, n) # propensity scores
#' M_set <- sample(1:n, size = n/2) # main set
#'
#' ## specify learner
#' learner <- "mlr3::lrn('ranger', num.trees = 10)"
#'
#' ## run single GenericML iteration
#' GenericML_single(Z, D, Y, learner, propensity_scores, M_set)
#' }
#'
#' @export
GenericML_single <- function(Z, D, Y,
                             learner,
                             propensity_scores,
                             M_set,
                             A_set                = setdiff(1:length(Y), M_set),
                             Z_CLAN               = NULL,
                             HT                   = FALSE,
                             quantile_cutoffs     = c(0.25, 0.5, 0.75),
                             X1_BLP               = setup_X1(),
                             X1_GATES             = setup_X1(),
                             diff_GATES           = setup_diff(),
                             diff_CLAN            = setup_diff(),
                             vcov_BLP             = setup_vcov(),
                             vcov_GATES           = setup_vcov(),
                             equal_variances_CLAN = FALSE,
                             significance_level   = 0.05,
                             min_variation        = 1e-05){

  # input checks
  InputChecks_D(D)
  InputChecks_Y(Y)
  InputChecks_Z(Z)
  InputChecks_Z_CLAN(Z_CLAN)
  InputChecks_equal.length3(D, Y, Z)
  N <- length(Y)
  stopifnot(is.character(learner))
  stopifnot(length(learner) == 1)
  InputChecks_index_set(A_set, N)
  InputChecks_index_set(M_set, N)
  InputChecks_X1(X1_BLP, N)
  InputChecks_X1(X1_GATES, N)
  InputChecks_vcov.control(vcov_BLP)
  InputChecks_vcov.control(vcov_GATES)
  InputChecks_diff(diff_GATES, K = length(quantile_cutoffs) + 1)
  InputChecks_diff(diff_CLAN, K = length(quantile_cutoffs) + 1)
  stopifnot(is.numeric(quantile_cutoffs))
  stopifnot(0 < min(quantile_cutoffs) & max(quantile_cutoffs) < 1)
  stopifnot(is.logical(equal_variances_CLAN))
  stopifnot(is.logical(HT))
  stopifnot(is.numeric(significance_level) & length(significance_level) == 1)
  stopifnot(0.0 < significance_level & significance_level < 0.5)
  stopifnot(is.numeric(min_variation) & min_variation > 0)
  stopifnot(is.character(learner))
  stopifnot(length(learner) == 1)
  stopifnot(is.numeric(propensity_scores))
  InputChecks_equal.length2(Y, propensity_scores)
  InputChecks_propensity_scores(propensity_scores)



  # if no input provided, set Z_CLAN equal to Z
  if(is.null(Z_CLAN)) Z_CLAN <- Z
  InputChecks_equal.length2(Z, Z_CLAN)

  # set variable names for CLAN
  if(is.null(colnames(Z_CLAN))) colnames(Z_CLAN) <- paste0("V", 1:ncol(Z_CLAN))
  if(any(colnames(Z_CLAN) == "")){

    idx <- which(colnames(Z_CLAN) == "")
    colnames(Z_CLAN)[idx] <- paste0("V", idx)

  } # IF

  # render the learner an mlr3 environment
  learner <- get.learner_regr(make.mlr3.environment(learner, regr = TRUE))


  # call the main function
  GenericML_single_NoChecks(Z = Z, D = D, Y = Y,
                            propensity_scores = propensity_scores,
                            learner = learner,
                            M_set = M_set, A_set = A_set,
                            Z_CLAN                     = Z_CLAN,
                            X1_BLP                     = X1_BLP,
                            X1_GATES                   = X1_GATES,
                            HT                         = HT,
                            vcov_BLP                   = setup_vcov_align(vcov_BLP),   # align for consistency
                            vcov_GATES                 = setup_vcov_align(vcov_GATES), # align for consistency
                            equal_variances_CLAN       = equal_variances_CLAN,
                            quantile_cutoffs           = quantile_cutoffs,
                            diff_GATES                 = diff_GATES,
                            diff_CLAN                  = diff_CLAN,
                            significance_level         = significance_level,
                            min_variation              = min_variation)

} # END FUN


# helper that skips the input checks
GenericML_single_NoChecks <-
  function(Z, D, Y,
           propensity_scores,
           learner,
           M_set, A_set,
           Z_CLAN                     = NULL,
           X1_BLP                     = setup_X1(),
           X1_GATES                   = setup_X1(),
           HT                         = FALSE,
           vcov_BLP                   = setup_vcov(),
           vcov_GATES                 = setup_vcov(),
           equal_variances_CLAN       = FALSE,
           quantile_cutoffs           = c(0.25, 0.5, 0.75),
           diff_GATES                 = setup_diff(),
           diff_CLAN                  = setup_diff(),
           significance_level         = 0.05,
           min_variation              = 1e-05){


    ### step 1a: learn proxy predictors by using the auxiliary set ----

    # get proxy baseline estimates
    proxy_BCA.obj <- proxy_BCA_NoChecks(
      Z = Z, D = D, Y = Y,
      A_set         = A_set,
      learner       = learner,
      min_variation = min_variation)

    # get estimates on main sample
    proxy_BCA_M     <- proxy_BCA.obj$estimates[M_set]

    # get the proxy estimates of the CATE
    proxy_CATE.obj <-
      proxy_CATE_NoChecks(
        Z = Z, D = D, Y = Y,
        A_set          = A_set,
        learner        = learner,
        proxy_BCA      = proxy_BCA.obj$estimates,
        min_variation  = min_variation)

    # get estimates on main sample
    proxy_CATE_M <- proxy_CATE.obj$estimates$CATE[M_set]

    # restrict the X1_controls to M_set
    X1_BLP_M   <-
      setup_X1_NoChecks(
       funs_Z        = X1_BLP$funs_Z,
       covariates    = X1_BLP$covariates[M_set,,drop = FALSE],
       fixed_effects = X1_BLP$fixed_effects[M_set])

    X1_GATES_M <-
      setup_X1_NoChecks(
        funs_Z        = X1_GATES$funs_Z,
        covariates    = X1_GATES$covariates[M_set,,drop = FALSE],
        fixed_effects = X1_GATES$fixed_effects[M_set])

    # restrict the vcov controls to M_set
    vcov_BLP_M   <- setup_vcov_subset(vcov_BLP, idx = M_set)
    vcov_GATES_M <- setup_vcov_subset(vcov_GATES, idx = M_set)


    ### step 1b: estimate BLP parameters ----
    blp.obj <- BLP_NoChecks(
      D                  = D[M_set],
      Y                  = Y[M_set],
      propensity_scores  = propensity_scores[M_set],
      proxy_BCA          = proxy_BCA_M,
      proxy_CATE         = proxy_CATE_M,
      HT                 = HT,
      X1_control         = X1_BLP_M,
      vcov_control       = vcov_BLP_M,
      significance_level = significance_level)


    ### step 1c: estimate GATES parameters by OLS ----
    # group the proxy estimators for the CATE in the main sample by quantiles
    membership_M <- quantile_group_NoChecks(proxy_CATE_M,
                                            cutoffs = quantile_cutoffs)

    # estimate GATES
    gates.obj <- GATES_NoChecks(
      D                   = D[M_set],
      Y                   = Y[M_set],
      propensity_scores   = propensity_scores[M_set],
      proxy_BCA           = proxy_BCA_M,
      proxy_CATE          = proxy_CATE_M,
      membership          = membership_M,
      HT                  = HT,
      X1_control          = X1_GATES_M,
      vcov_control        = vcov_GATES_M,
      diff                = diff_GATES,
      significance_level  = significance_level)


    ### step 1d: estimate CLAN parameters on the main sample ----
    clan.obj <- CLAN_NoChecks(
      Z_CLAN             = Z_CLAN[M_set,,drop = FALSE],
      membership         = membership_M,
      equal_variances    = equal_variances_CLAN,
      diff               = diff_CLAN,
      significance_level = significance_level)


    ### step 1e: get parameters over which we maximize to find the "best" ML method ----
    best.obj <- lambda_parameters_NoChecks(BLP        = blp.obj,
                                           GATES      = gates.obj,
                                           proxy_CATE = proxy_CATE_M,
                                           membership = membership_M)

    ### organize output in a list ----
    return(list(BLP            = blp.obj,
                GATES          = gates.obj,
                CLAN           = clan.obj,
                proxy_CATE     = proxy_CATE.obj,
                proxy_BCA      = proxy_BCA.obj,
                best           = best.obj
    ))

  } # FUN



#' Estimate the two lambda parameters
#'
#' Estimates the lambda parameters \eqn{\Lambda} and \eqn{\bar{\Lambda}} whose medians are used to find the best ML method.
#'
#' @param BLP An object of class \code{"\link{BLP}"}.
#' @param GATES An object of class \code{"\link{GATES}"}.
#' @param proxy_CATE Proxy estimates of the CATE.
#' @param membership A logical matrix that indicates the group membership of each observation in \code{Z_CLAN}. Needs to be of type \code{"\link{quantile_group}"}. Typically, the grouping is based on CATE estimates, which are for instance returned by \code{\link{proxy_CATE}()}.
#'
#' @return
#' A list containing the estimates of \eqn{\Lambda} and \eqn{\bar{\Lambda}}, denoted \code{lambda} and \code{lambda.bar}, respectively.
#'
#' @references
#' Chernozhukov V., Demirer M., Duflo E., Fernández-Val I. (2020). \dQuote{Generic Machine Learning Inference on Heterogenous Treatment Effects in Randomized Experiments.} \emph{arXiv preprint arXiv:1712.04802}. URL: \url{https://arxiv.org/abs/1712.04802}.
#'
#' @examples
#' ## generate data
#' set.seed(1)
#' n  <- 200                                # number of observations
#' p  <- 5                                  # number of covariates
#' D  <- rbinom(n, 1, 0.5)                  # random treatment assignment
#' Y  <- runif(n)                           # outcome variable
#' propensity_scores <- rep(0.5, n)         # propensity scores
#' proxy_BCA         <- runif(n)            # proxy BCA estimates
#' proxy_CATE        <- runif(n)            # proxy CATE estimates
#' membership <- quantile_group(proxy_CATE) # group membership
#'
#' ## perform BLP
#' BLP <- BLP(Y, D, propensity_scores, proxy_BCA, proxy_CATE)
#'
#' ## perform GATES
#' GATES <- GATES(Y, D, propensity_scores, proxy_BCA, proxy_CATE, membership)
#'
#' ## get estimates of the lambda parameters
#' lambda_parameters(BLP, GATES, proxy_CATE, membership)
#'
#' @export
lambda_parameters <- function(BLP,
                              GATES,
                              proxy_CATE,
                              membership){

  if(!inherits(x = BLP, what = "BLP", which = FALSE)){
    stop("The BLP object must be an instance of BLP()")
  }
  if(!inherits(x = GATES, what = "GATES", which = FALSE)){
    stop("The GATES object must be an instance of GATES()")
  }
  InputChecks_group.membership(membership)
  stopifnot(is.numeric(proxy_CATE))
  InputChecks_equal.length2(proxy_CATE, membership)
  temp <- GATES$coefficients
  gates.coefs <- temp[startsWith(rownames(temp), "gamma."), "Estimate"]

  if(ncol(membership) != length(gates.coefs)){
    stop("The number of columns of 'membership' must be equal to the number of GATES gamma coefficients")
  }

  # return
  lambda_parameters_NoChecks(BLP = BLP,
                             GATES = GATES,
                             proxy_CATE = proxy_CATE,
                             membership = membership)

} # END FUN


# same as above, but w/o input checks
lambda_parameters_NoChecks <- function(BLP,
                                       GATES,
                                       proxy_CATE,
                                       membership){

  temp <- GATES$coefficients
  gates.coefs <- temp[startsWith(rownames(temp), "gamma."), "Estimate"]

  return(list(lambda = BLP$coefficients["beta.2", "Estimate"]^2 * stats::var(proxy_CATE),
              lambda.bar = as.numeric(colMeans(membership) %*%  gates.coefs^2)))

} # FUN
