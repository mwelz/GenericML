#' Performs BLP regression
#'
#' Performs the linear regression for the Best Linear Predictor (BLP) procedure.
#'
#' @param Y A numeric vector containing the response variable.
#' @param D A binary vector of treatment assignment. Value one denotes assignment to the treatment group and value zero assignment to the control group.
#' @param propensity_scores A numeric vector of propensity scores. We recommend to use the estimates of a \code{\link{propensity_score}} object.
#' @param proxy_BCA A numeric vector of proxy baseline conditional average (BCA) estimates. We recommend to use the estimates of a \code{\link{proxy_BCA}} object.
#' @param proxy_CATE A numeric vector of proxy conditional average treatment effect (CATE) estimates. We recommend to use the estimates of a \code{\link{proxy_CATE}} object.
#' @param HT Logical. If \code{TRUE}, a Horvitz-Thompson (HT) transformation is applied (BLP2 in the paper). Default is \code{FALSE}.
#' @param X1_control Specifies the design matrix \eqn{X_1} in the regression. Must be an instance of \code{\link{setup_X1}}. See the documentation of \code{\link{setup_X1}} for details.
#' @param vcov_control Specifies the covariance matrix estimator. Must be an instance of \code{\link{setup_vcov}}. See the documentation of \code{\link{setup_vcov}} for details.
#' @param significance_level Significance level. Default is 0.05.
#'
#' @return
#' An object of class \code{BLP}, consisting of the following components:
#' \describe{
#'   \item{\code{generic_targets}}{A matrix of the inferential results on the BLP generic targets.}
#'   \item{\code{coefficients}}{An object of class \code{\link[lmtest]{coeftest}}, contains the coefficients of the BLP regression.}
#'   \item{\code{lm}}{An object of class \code{\link[stats]{lm}} used to fit the linear regression model.}
#'   }
#'
#' @references
#' Chernozhukov V., Demirer M., Duflo E., Fernández-Val I. (2020). \dQuote{Generic Machine Learning Inference on Heterogenous Treatment Effects in Randomized Experiments.} \emph{arXiv preprint arXiv:1712.04802}. URL: \url{https://arxiv.org/abs/1712.04802}.
#'
#' @seealso
#' \code{\link{setup_X1}},
#' \code{\link{setup_diff}},
#' \code{\link{setup_vcov}},
#' \code{\link{propensity_score}},
#' \code{\link{proxy_BCA}},
#' \code{\link{proxy_CATE}}
#'
#' @examples
#' ## generate data
#' set.seed(1)
#' n  <- 150                        # number of observations
#' p  <- 5                          # number of covariates
#' D  <- rbinom(n, 1, 0.5)          # random treatment assignment
#' Y  <- runif(n)                   # outcome variable
#' propensity_scores <- rep(0.5, n) # propensity scores
#' proxy_BCA         <- runif(n)    # proxy BCA estimates
#' proxy_CATE        <- runif(n)    # proxy CATE estimates
#'
#' ## perform BLP
#' BLP(Y, D, propensity_scores, proxy_BCA, proxy_CATE)
#'
#' @export
BLP <- function(Y, D,
                propensity_scores,
                proxy_BCA,
                proxy_CATE,
                HT                 = FALSE,
                X1_control         = setup_X1(),
                vcov_control       = setup_vcov(),
                significance_level = 0.05){


  # input checks
  InputChecks_Y(Y)
  InputChecks_D(D)
  InputChecks_equal.length2(Y, D)
  InputChecks_equal.length2(proxy_BCA, proxy_CATE)
  InputChecks_equal.length2(proxy_BCA, Y)
  InputChecks_vcov.control(vcov_control)
  InputChecks_X1(X1_control, length(Y))
  stopifnot(is.numeric(proxy_BCA))
  stopifnot(is.numeric(proxy_CATE))
  stopifnot(is.logical(HT))
  stopifnot(is.numeric(significance_level) & length(significance_level) == 1)
  stopifnot(0.0 < significance_level & significance_level < 0.5)
  stopifnot(is.numeric(propensity_scores))
  InputChecks_equal.length2(Y, propensity_scores)
  InputChecks_propensity_scores(propensity_scores)

  # fit model according to strategy 1 or 2 in the paper
  BLP_NoChecks(D = D, Y = Y,
               propensity_scores  = propensity_scores,
               proxy_BCA     = proxy_BCA,
               proxy_CATE         = proxy_CATE,
               X1_control         = X1_control,
               vcov_control       = vcov_control,
               significance_level = significance_level)

} # FUN


# helper function to perform BLP w/o input checks
BLP_NoChecks <- function(D, Y,
                         propensity_scores,
                         proxy_BCA,
                         proxy_CATE,
                         HT                 = FALSE,
                         X1_control         = setup_X1(),
                         vcov_control       = setup_vcov(),
                         significance_level = 0.05){

  # fit model according to strategy 1 or 2 in the paper
  do.call(what = get(ifelse(HT, "BLP.HT", "BLP.classic")),
          args = list(D = D, Y = Y,
                      propensity_scores  = propensity_scores,
                      proxy_BCA     = proxy_BCA,
                      proxy_CATE         = proxy_CATE,
                      X1_control         = X1_control,
                      vcov_control       = vcov_control,
                      significance_level = significance_level))

} # FUN


# helper function for case when there is no HT transformation used. Wrapped by function "BLP"
BLP.classic <- function(D, Y, propensity_scores,
                        proxy_BCA, proxy_CATE,
                        X1_control         = setup_X1(),
                        vcov_control       = setup_vcov(),
                        significance_level = 0.05){

  # prepare weights
  weights <- 1 / (propensity_scores * (1 - propensity_scores))

  # prepare covariate matrix X
  X <- data.frame(get.df.from.X1_control(functions.of.Z_mat = cbind(S = proxy_CATE,
                                                                    B = proxy_BCA,
                                                                    p = propensity_scores),
                                           X1_control = X1_control),
                  beta.1 = D - propensity_scores,
                  beta.2 = (D - propensity_scores) * (proxy_CATE - mean(proxy_CATE)))

  # fit weighted linear regression by OLS
  blp.obj <- stats::lm(Y ~., data = data.frame(Y, X), weights = weights)

  # get estimate of covariance matrix of the error terms
  vcov. <- get.vcov(x              = blp.obj,
                    vcov_control   = vcov_control)

  # extract coefficients
  coefficients <- lmtest::coeftest(blp.obj, vcov. = vcov.)

  # return
  return(structure(
    list(generic_targets = generic_targets_BLP(coeftest.object = coefficients,
                                               significance_level = significance_level),
         coefficients = coefficients,
         lm = blp.obj),
    class = "BLP"))

} # END FUN



# helper function for case when there is a HT transformation used. Wrapped by function "BLP"
BLP.HT <- function(D, Y, propensity_scores,
                   proxy_BCA, proxy_CATE,
                   X1_control         = setup_X1(),
                   vcov_control       = setup_vcov(),
                   significance_level = 0.05){

  # HT transformation
  H <- (D - propensity_scores) / (propensity_scores * (1 - propensity_scores))

  # prepare matrix X1
  X1. <- get.df.from.X1_control(functions.of.Z_mat = cbind(S = proxy_CATE,
                                                           B = proxy_BCA,
                                                           p = propensity_scores),
                                 X1_control = X1_control)


  # construct the matrix X1H (the fixed effects are not multiplied by H, if applicable)
  if(is.null(X1_control$fixed_effects)){

    # matrix X_1 * H
    X1H           <- X1. * H
    colnames(X1H) <- paste0(colnames(X1.), ".H")

  } else{

    fixed.effects     <- X1.$fixed.effects  # retain the fixed effects
    X1.               <- X1.[,!colnames(X1.) %in% "fixed.effects", drop = FALSE]
    X1H               <- X1. * H
    colnames(X1H)     <- paste0(colnames(X1.), ".H")
    X1H$fixed.effects <- fixed.effects # append the fixed effects

  } # IF

  # prepare covariate matrix X
  X <- data.frame(X1H,
                  beta.2 = proxy_CATE - mean(proxy_CATE))

  # fit linear regression by OLS (intercept is beta.1)
  blp.obj <- stats::lm(YH ~., data = data.frame(YH = Y*H, X))
  names(blp.obj$coefficients) <- c("beta.1", names(blp.obj$coefficients)[-1])

  # get estimate of covariance matrix of the error terms
  vcov. <- get.vcov(x              = blp.obj,
                    vcov_control   = vcov_control)

  # extract coefficients
  coefficients <- lmtest::coeftest(blp.obj, vcov. = vcov.)

  # return
  return(structure(
    list(generic_targets = generic_targets_BLP(coeftest.object = coefficients,
                                               significance_level = significance_level),
         coefficients = coefficients,
         lm = blp.obj),
    class = "BLP"))

} # END FUN


# helper function to calculate the generic targets of BLP
generic_targets_BLP <- function(coeftest.object, significance_level = 0.05){

  # helper that computes generic targets beta.1, beta.2 in BLP
  coefficients.temp <- coeftest.object[c("beta.1", "beta.2"), 1:3]
  colnames(coefficients.temp) <- c("Estimate", "Std. Error", "z value")
  p.right <- stats::pnorm(coefficients.temp[,"z value"], lower.tail = FALSE) # right p-value: Pr(Z>z)
  p.left  <- stats::pnorm(coefficients.temp[,"z value"], lower.tail = TRUE)  # left p-value: Pr(Z<z)
  ci.lo   <- coefficients.temp[,"Estimate"] - stats::qnorm(1-significance_level/2) * coefficients.temp[,"Std. Error"]
  ci.up   <- coefficients.temp[,"Estimate"] + stats::qnorm(1-significance_level/2) * coefficients.temp[,"Std. Error"]
  generic_targets <- cbind(coefficients.temp, ci.lo, ci.up, p.left, p.right)
  colnames(generic_targets) <- c("Estimate", "Std. Error", "z value",
                                 "CB lower", "CB upper", "Pr(<z)", "Pr(>z)")
  generic_targets <- generic_targets[,c("Estimate", "CB lower", "CB upper",
                                        "Std. Error", "z value", "Pr(<z)", "Pr(>z)")]

  return(generic_targets)

} # FUN
