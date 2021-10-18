#' Estimates the BLP parameters based on the main sample M.
#'
#' @param D a binary vector of treatment status of length _|M|_
#' @param Y a vector of responses of length _|M|_
#' @param propensity_scores a vector of propensity scores of length _|M|_
#' @param proxy.baseline a vector of proxy baseline estimates of length _M_
#' @param proxy.cate a vector of proxy CATE estimates of length _M_
#' @param HT logical. If TRUE, a HT transformation is applied (BLP2 in the paper). Default is FALSE.
#' @param X1_control Specifies the design matrix \eqn{X_1} in the regression. See the documentation of \code{\link{setup_X1}} for details.
#' @param vcov_control Specifies the covariance matrix estimator. See the documentation of \code{\link{setup_vcov}} for details.
#' @param significance_level significance level for construction of confidence intervals
#' @return BLP coefficients with inference statements
#'
#' @export
BLP <- function(D, Y,
                propensity_scores,
                proxy.baseline,
                proxy.cate,
                HT  = FALSE,
                X1_control       = setup_X1(),
                vcov_control       = setup_vcov(),
                significance_level = 0.05){


  # input checks
  InputChecks_Y(Y)
  InputChecks_D(D)
  InputChecks_equal.length2(Y, D)
  InputChecks_equal.length2(proxy.baseline, proxy.cate)
  InputChecks_equal.length2(proxy.baseline, Y)
  InputChecks_vcov.control(vcov_control)
  InputChecks_X1(X1_control)

  # fit model according to strategy 1 or 2 in the paper
  BLP_NoChecks(D = D, Y = Y,
               propensity_scores  = propensity_scores,
               proxy.baseline     = proxy.baseline,
               proxy.cate         = proxy.cate,
               X1_control       = X1_control,
               vcov_control       = vcov_control,
               significance_level = significance_level)

} # FUN


# helper function to perform BLP w/o input checks
BLP_NoChecks <- function(D, Y,
                         propensity_scores,
                         proxy.baseline,
                         proxy.cate,
                         HT  = FALSE,
                         X1_control       = setup_X1(),
                         vcov_control       = setup_vcov(),
                         significance_level = 0.05){

  # fit model according to strategy 1 or 2 in the paper
  do.call(what = get(ifelse(HT, "BLP.HT", "BLP.classic")),
          args = list(D = D, Y = Y,
                      propensity_scores  = propensity_scores,
                      proxy.baseline     = proxy.baseline,
                      proxy.cate         = proxy.cate,
                      X1_control       = X1_control,
                      vcov_control       = vcov_control,
                      significance_level = significance_level))

} # FUN


# helper function for case when there is no HT transformation used. Wrapped by function "BLP"
BLP.classic <- function(D, Y, propensity_scores,
                        proxy.baseline, proxy.cate,
                        X1_control       = setup_X1(),
                        vcov_control       = setup_vcov(),
                        significance_level = 0.05){

  # prepare weights
  weights <- 1 / (propensity_scores * (1 - propensity_scores))

  # prepare covariate matrix X
  X <- data.frame(get.df.from.X1_control(functions.of.Z_mat = cbind(S = proxy.cate,
                                                                      B = proxy.baseline,
                                                                      p = propensity_scores),
                                           X1_control = X1_control),
                  beta.1 = D - propensity_scores,
                  beta.2 = (D - propensity_scores) * (proxy.cate - mean(proxy.cate)))

  # fit weighted linear regression by OLS
  blp.obj <- stats::lm(Y ~., data = data.frame(Y, X), weights = weights)

  # get estimate of covariance matrix of the error terms
  vcov. <- get.vcov(x              = blp.obj,
                    vcov_control   = vcov_control)

  # extract coefficients
  coefficients <- lmtest::coeftest(blp.obj, vcov. = vcov.)

  # return
  return(structure(
    list(lm.obj = blp.obj,
              blp.coefficients = blp.obj$coefficients[c("beta.1", "beta.2")],
              generic.targets = generic.targets_BLP(coefficients, significance_level = significance_level),
              coefficients = coefficients), class = "BLP"))

} # END FUN



# helper function for case when there is a HT transformation used. Wrapped by function "BLP"
BLP.HT <- function(D, Y, propensity_scores,
                   proxy.baseline, proxy.cate,
                   X1_control       = setup_X1(),
                   vcov_control       = setup_vcov(),
                   significance_level = 0.05){

  # HT transformation
  H <- (D - propensity_scores) / (propensity_scores * (1 - propensity_scores))

  # prepare matrix X1
  X1. <- get.df.from.X1_control(functions.of.Z_mat = cbind(S = proxy.cate,
                                                             B = proxy.baseline,
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
                  beta.2 = proxy.cate - mean(proxy.cate))

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
    list(lm.obj = blp.obj,
              blp.coefficients = blp.obj$coefficients[c("beta.1", "beta.2")],
              generic.targets = generic.targets_BLP(coefficients, significance_level = significance_level),
              coefficients = coefficients), class = "BLP"))

} # END FUN


# helper function to calculate the generic targets of BLP
generic.targets_BLP <- function(coeftest.object, significance_level = 0.05){

  # helper that computes generic targets beta.1, beta.2 in BLP
  coefficients.temp <- coeftest.object[c("beta.1", "beta.2"), 1:3]
  colnames(coefficients.temp) <- c("Estimate", "Std. Error", "z value")
  p.right <- stats::pnorm(coefficients.temp[,"z value"], lower.tail = FALSE) # right p-value: Pr(Z>z)
  p.left  <- stats::pnorm(coefficients.temp[,"z value"], lower.tail = TRUE)  # left p-value: Pr(Z<z)
  ci.lo   <- coefficients.temp[,"Estimate"] - stats::qnorm(1-significance_level/2) * coefficients.temp[,"Std. Error"]
  ci.up   <- coefficients.temp[,"Estimate"] + stats::qnorm(1-significance_level/2) * coefficients.temp[,"Std. Error"]
  generic.targets <- cbind(coefficients.temp, ci.lo, ci.up, p.left, p.right)
  colnames(generic.targets) <- c("Estimate", "Std. Error", "z value",
                                 "CB lower", "CB upper", "Pr(<z)", "Pr(>z)")
  generic.targets <- generic.targets[,c("Estimate", "CB lower", "CB upper",
                                        "Std. Error", "z value", "Pr(<z)", "Pr(>z)")]

  return(generic.targets)

} # FUN
