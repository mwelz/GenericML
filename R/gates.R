#' Estimates the GATES parameters based on the main sample M.
#'
#' @param D a binary vector of treatment status of length _|M|_
#' @param Y a vector of responses of length _|M|_
#' @param propensity.scores a vector of propensity scores of length _|M|_
#' @param proxy.baseline a vector of proxy baseline estimates of length _M_
#' @param proxy.cate a vector of proxy CATE estimates of length _M_
#' @param group.membership.main.sample a logical matrix with _M_ rows that indicate
#' the group memberships (such a matrix is returned by the function quantile_group())
#' @param HT logical. If TRUE, a HT transformation is applied (GATES2 in the paper). Default is FALSE.
#' @param X1_control Specifies the design matrix \eqn{X_1} in the regression. See the documentation of \code{\link{setup_X1}} for details.
#' @param vcov_control Specifies the covariance matrix estimator. See the documentation of \code{\link{setup_vcov}} for details.
#' @param diff Controls the generic targets of GATES See the documentation of \code{\link{setup_diff}}.
#' @param significance_level significance level for construction of confidence intervals
#'
#' @export
GATES <- function(D, Y,
                  propensity.scores,
                  proxy.baseline,
                  proxy.cate,
                  group.membership.main.sample,
                  HT                  = FALSE,
                  X1_control        = setup_X1(),
                  vcov_control        = setup_vcov(),
                  diff = setup_diff(),
                  significance_level  = 0.05){

  # input check
  InputChecks_D(D)
  InputChecks_Y(Y)
  InputChecks_equal.length2(D, Y)
  InputChecks_equal.length3(propensity.scores, proxy.baseline, proxy.cate)
  InputChecks_equal.length2(Y, propensity.scores)
  InputChecks_X1(X1_control)
  InputChecks_vcov.control(vcov_control)
  InputChecks_diff(diff, K = ncol(group.membership.main.sample))
  InputChecks_group.membership(group.membership.main.sample)

  # fit model according to strategy 1 or 2 in the paper
  GATES_NoChecks(D = D, Y = Y,
                 propensity.scores   = propensity.scores,
                 proxy.baseline      = proxy.baseline,
                 proxy.cate          = proxy.cate,
                 group.membership.main.sample = group.membership.main.sample,
                 X1_control        = X1_control,
                 vcov_control        = vcov_control,
                 diff = diff,
                 significance_level  = significance_level)

} # FUN


# helper function that skips the input checks
GATES_NoChecks <- function(D, Y,
                           propensity.scores,
                           proxy.baseline,
                           proxy.cate,
                           group.membership.main.sample,
                           HT                  = FALSE,
                           X1_control        = setup_X1(),
                           vcov_control        = setup_vcov(),
                           diff                = setup_diff(),
                           significance_level  = 0.05){

  # fit model according to strategy 1 or 2 in the paper
  do.call(what = get(ifelse(HT, "GATES.HT", "GATES.classic")),
          args = list(D = D, Y = Y,
                      propensity.scores   = propensity.scores,
                      proxy.baseline      = proxy.baseline,
                      proxy.cate          = proxy.cate,
                      group.membership.main.sample = group.membership.main.sample,
                      X1_control        = X1_control,
                      vcov_control        = vcov_control,
                      diff                = diff,
                      significance_level  = significance_level))

} # FUN


# helper function for case when there is no HT transformation used. Wrapped by function "GATES"
GATES.classic <- function(D, Y,
                          propensity.scores,
                          proxy.baseline, proxy.cate,
                          group.membership.main.sample,
                          X1_control       = setup_X1(),
                          vcov_control       = setup_vcov(),
                          diff               = setup_diff(),
                          significance_level = 0.05){

  # make the group membership a binary matrix
  groups <- 1 * group.membership.main.sample

  # number of groups
  K <- ncol(groups)

  # prepare weights
  weights <- 1 / (propensity.scores * (1 - propensity.scores))

  # prepare matrix X1
  X1     <- get.df.from.X1_control(functions.of.Z_mat = cbind(S = proxy.cate,
                                                                B = proxy.baseline,
                                                                p = propensity.scores),
                                     X1_control = X1_control)

  # prepare covariate matrix
  X <- data.frame(X1,
                  (D - propensity.scores) * groups)
  colnames(X) <- c(colnames(X1), paste0("gamma.", 1:K))

  # fit weighted linear regression by OLS
  gates.obj <- stats::lm(Y ~., data = data.frame(Y, X), weights = weights)

  # get estimate of covariance matrix of the error terms
  vcov. <- get.vcov(x              = gates.obj,
                    vcov_control   = vcov_control)

  # extract the relevant coefficients
  coefficients                 <- lmtest::coeftest(gates.obj, vcov. = vcov.)
  gates.coefficients           <- coefficients[paste0("gamma.", 1:K), 1]
  gates.coefficients.quantiles <- colnames(groups)
  names(gates.coefficients.quantiles) <- paste0("gamma.", 1:K)

  # return
  return(structure(
    list(lm.obj = gates.obj,
              gates.coefficients = gates.coefficients,
              gates.coefficients.quantiles = gates.coefficients.quantiles,
              generic.targets = generic.targets_GATES(coeftest.object = coefficients,
                                                      K = K,
                                                      vcov = vcov.,
                                                      significance_level = significance_level,
                                                      diff = diff),
              coefficients = coefficients), class = "GATES"))

} # END FUN


# helper function for case when there is a HT transformation used. Wrapped by function "GATES"
GATES.HT <- function(D, Y,
                     propensity.scores,
                     proxy.baseline, proxy.cate,
                     group.membership.main.sample,
                     X1_control       = setup_X1(),
                     vcov_control       = setup_vcov(),
                     diff               = setup_diff(),
                     significance_level = 0.05){

  # make the group membership a binary matrix
  groups <- 1 * group.membership.main.sample

  # number of groups
  K <- ncol(groups)

  # HT transformation
  H <- (D - propensity.scores) / (propensity.scores * (1 - propensity.scores))

  # prepare matrix X1
  X1 <- get.df.from.X1_control(functions.of.Z_mat = cbind(S = proxy.cate,
                                                            B = proxy.baseline,
                                                            p = propensity.scores),
                                 X1_control = X1_control)

  # construct the matrix X1H (the fixed effects are not multiplied by H, if applicable)
  if(is.null(X1_control$fixed_effects)){

    # matrix X_1 * H
    X1H           <- X1 * H
    colnames(X1H) <- paste0(colnames(X1), ".H")

  } else{

    fixed.effects     <- X1$fixed.effects  # retain the fixed effects
    X1                <- X1[,!colnames(X1) %in% "fixed.effects", drop = FALSE]
    X1H               <- X1 * H
    colnames(X1H)     <- paste0(colnames(X1), ".H")

    # perform one-hot encoding on fixed effect categories and drop first category to avoid dummy trap
    X1H <- data.frame(X1H,
                      as.matrix(stats::model.matrix(~ fixed.effects + 0, data = fixed.effects))[,-1])

  } # IF


  # prepare covariate matrix
  X <- data.frame(X1H,  groups)
  colnames(X) <- c(colnames(X1H), paste0("gamma.", 1:K))

  # fit linear regression by OLS (no intercept!)
  gates.obj <- stats::lm(formula =  stats::as.formula(paste0("YH ~ ", paste0(colnames(X), collapse = " + "), " + 0")),
                data = data.frame(YH = Y*H, X))

  # get estimate of covariance matrix of the error terms
  vcov. <- get.vcov(x              = gates.obj,
                    vcov_control   = vcov_control)

  # extract the relevant coefficients
  coefficients                 <- lmtest::coeftest(gates.obj, vcov. = vcov.)
  gates.coefficients           <- coefficients[paste0("gamma.", 1:K), 1]
  gates.coefficients.quantiles <- colnames(groups)
  names(gates.coefficients.quantiles) <- paste0("gamma.", 1:K)

  # return
  return(structure(
    list(lm.obj = gates.obj,
              gates.coefficients = gates.coefficients,
              gates.coefficients.quantiles = gates.coefficients.quantiles,
              generic.targets = generic.targets_GATES(coeftest.object = coefficients,
                                                      K = K,
                                                      vcov = vcov.,
                                                      significance_level = significance_level,
                                                      diff = diff),
              coefficients = coefficients), class = "GATES"))

} # END FUN


# helper function to calculate the generic targets of BLP (vcov is estimate, not list)
generic.targets_GATES <- function(coeftest.object, K, vcov,
                                  significance_level = 0.05,
                                  diff = diff){

  # extract controls
  subtract_from  <- diff$subtract_from
  subtracted <- diff$subtracted

  # extract coefficients
  coefficients.temp <- coeftest.object[paste0("gamma.", 1:K), 1:3]
  colnames(coefficients.temp) <- c("Estimate", "Std. Error", "z value")

  # get quantile
  z          <- stats::qnorm(1-significance_level/2)

  # compute relevant statistics
  p.right <- stats::pnorm(coefficients.temp[,"z value"], lower.tail = FALSE) # right p-value: Pr(Z>z)
  p.left  <- stats::pnorm(coefficients.temp[,"z value"], lower.tail = TRUE)  # left p-value: Pr(Z<z)
  ci.lo   <- coefficients.temp[,"Estimate"] - z * coefficients.temp[,"Std. Error"]
  ci.up   <- coefficients.temp[,"Estimate"] + z * coefficients.temp[,"Std. Error"]

  # prepare generic targets of the gammas
  generic.targets <- cbind(coefficients.temp, ci.lo, ci.up, p.left, p.right)
  colnames(generic.targets) <- c("Estimate", "Std. Error", "z value",
                                 "CB lower", "CB upper", "Pr(<z)", "Pr(>z)")
  generic.targets <- generic.targets[,c("Estimate", "CB lower", "CB upper",
                                        "Std. Error", "z value", "Pr(<z)", "Pr(>z)")]

  # GATES differences: estimate mean and variance of a difference
  if(subtract_from == "least"){

    diff. <- coeftest.object["gamma.1", "Estimate"] -
      coeftest.object[paste0("gamma.", subtracted), "Estimate"]

    diff.se <- as.numeric(
      sqrt(vcov["gamma.1", "gamma.1"] +
                      diag(vcov[paste0("gamma.", subtracted),
                                paste0("gamma.", subtracted), drop = FALSE]) -
                      2 * vcov["gamma.1", paste0("gamma.", subtracted)]))

    nam <- paste0("gamma.1-gamma.", subtracted)

  } else{

    diff. <- coeftest.object[paste0("gamma.", K), "Estimate"] -
      coeftest.object[paste0("gamma.", subtracted), "Estimate"]

    diff.se <- as.numeric(
      sqrt(vcov[paste0("gamma.", K), paste0("gamma.", K)] +
                      diag(vcov[paste0("gamma.", subtracted),
                                paste0("gamma.", subtracted), drop = FALSE]) -
                      2 * vcov[paste0("gamma.", K), paste0("gamma.", subtracted)]))

    nam <- paste0("gamma.", K, "-gamma.", subtracted)

  } # IF

  ci.lo   <- diff. - z * diff.se
  ci.up   <- diff. + z * diff.se
  zstat   <- diff. / diff.se
  p.right <- stats::pnorm(zstat, lower.tail = FALSE) # right p-value: Pr(Z>z)
  p.left  <- stats::pnorm(zstat, lower.tail = TRUE)  # left p-value: Pr(Z<z)

  diff.mat <- cbind(diff., ci.lo, ci.up, diff.se, zstat, p.left, p.right)
  rownames(diff.mat) <- nam

  # return final matrix
  rbind(generic.targets, diff.mat)

} # FUN
