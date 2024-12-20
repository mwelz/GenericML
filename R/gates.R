#' Performs GATES regression
#'
#' Performs the linear regression for the Group Average Treatments Effects (GATES) procedure.
#'
#' @param Y A numeric vector containing the response variable.
#' @param D A binary vector of treatment assignment. Value one denotes assignment to the treatment group and value zero assignment to the control group.
#' @param propensity_scores A numeric vector of propensity scores. We recommend to use the estimates of a \code{"\link{propensity_score}"} object.
#' @param proxy_BCA A numeric vector of proxy baseline conditional average (BCA) estimates. We recommend to use the estimates of a \code{"\link{proxy_BCA}"} object.
#' @param proxy_CATE A numeric vector of proxy conditional average treatment effect (CATE) estimates. We recommend to use the estimates of a \code{"\link{proxy_CATE}"} object.
#' @param membership A logical matrix that indicates the group membership of each observation in \code{Z_CLAN}. Needs to be of type \code{"\link{quantile_group}"}. Typically, the grouping is based on CATE estimates, which are for instance returned by \code{\link{proxy_CATE}()}.
#' @param HT Logical. If \code{TRUE}, a Horvitz-Thompson (HT) transformation is applied (GATES2 in the paper). Default is \code{FALSE}.
#' @param X1_control Specifies the design matrix \eqn{X_1} in the regression. Must be an object of class  \code{"\link{setup_X1}"}. See the documentation of \code{\link{setup_X1}()} for details.
#' @param vcov_control Specifies the covariance matrix estimator. Must be an object of class \code{"\link{setup_vcov}"}. See the documentation of \code{\link{setup_vcov}()} for details.
#' @param diff Specifies the generic targets of CLAN. Must be an object of class \code{"\link{setup_diff}"}. See the documentation of \code{\link{setup_diff}()} for details.
#' @param monotonize Logical. Should GATES point estimates and confidence bounds be rearranged to be monotonically increasing following the monotonization method of Chernozhukov et al. (2009, Biometrika)? Default is \code{TRUE}.
#' @param external_weights Optional vector of external numeric weights for weighted regression (in addition to the standard weights used when \code{HT = FALSE}).
#' @param significance_level Significance level. Default is 0.05.
#'
#' @return
#' An object of class \code{"GATES"}, consisting of the following components:
#' \describe{
#'   \item{\code{generic_targets}}{A matrix of the inferential results on the GATES generic targets.}
#'   \item{\code{coefficients}}{An object of class \code{"\link[lmtest]{coeftest}"}, contains the coefficients of the GATES regression.}
#'   \item{\code{lm}}{An object of class \code{"\link[stats]{lm}"} used to fit the linear regression model.}
#'   }
#'
#' @references
#' Chernozhukov V., Demirer M., Duflo E., Fernández-Val I. (2020). \dQuote{Generic Machine Learning Inference on Heterogenous Treatment Effects in Randomized Experiments.} \emph{arXiv preprint arXiv:1712.04802}. URL: \url{https://arxiv.org/abs/1712.04802}.
#' Chernozhukov V., Fernández-Val I., Galichon, A. (2009). \dQuote{Improving Point and Interval Estimators of Monotone Functions by Rearrangement.} \emph{Biometrika}, \bold{96}(3), 559--575. \doi{10.1093/biomet/asp030}.
#'
#' @examples
#' ## generate data
#' set.seed(1)
#' n  <- 150                                # number of observations
#' p  <- 5                                  # number of covariates
#' D  <- rbinom(n, 1, 0.5)                  # random treatment assignment
#' Y  <- runif(n)                           # outcome variable
#' propensity_scores <- rep(0.5, n)         # propensity scores
#' proxy_BCA         <- runif(n)            # proxy BCA estimates
#' proxy_CATE        <- runif(n)            # proxy CATE estimates
#' membership <- quantile_group(proxy_CATE) # group membership
#'
#' ## perform GATES
#' GATES(Y, D, propensity_scores, proxy_BCA, proxy_CATE, membership)
#'
#' @seealso
#' \code{\link{setup_X1}()},
#' \code{\link{setup_diff}()},
#' \code{\link{setup_vcov}()},
#' \code{\link{propensity_score}()},
#' \code{\link{proxy_BCA}()},
#' \code{\link{proxy_CATE}()}
#'
#' @export
GATES <- function(Y, D,
                  propensity_scores,
                  proxy_BCA,
                  proxy_CATE,
                  membership,
                  HT                 = FALSE,
                  X1_control         = setup_X1(),
                  vcov_control       = setup_vcov(),
                  diff               = setup_diff(),
                  monotonize         = TRUE,
                  external_weights   = NULL,
                  significance_level = 0.05){

  # input check
  InputChecks_D(D)
  InputChecks_Y(Y)
  InputChecks_equal.length2(D, Y)
  InputChecks_equal.length3(propensity_scores, proxy_BCA, proxy_CATE)
  stopifnot(is.numeric(propensity_scores))
  InputChecks_equal.length2(Y, propensity_scores)
  InputChecks_propensity_scores(propensity_scores)
  InputChecks_X1(X1_control, length(Y))
  InputChecks_vcov.control(vcov_control)
  InputChecks_diff(diff, K = ncol(membership))
  InputChecks_group.membership(membership)
  stopifnot(is.numeric(proxy_BCA))
  stopifnot(is.numeric(proxy_CATE))
  stopifnot(is.logical(HT))
  stopifnot(is.logical(monotonize))
  stopifnot(is.numeric(significance_level) & length(significance_level) == 1)
  stopifnot(0.0 < significance_level & significance_level < 0.5)
  InputChecks_external_weights(external_weights, length(Y))

  # fit model according to strategy 1 or 2 in the paper
  GATES_NoChecks(D = D, Y = Y,
                 propensity_scores   = propensity_scores,
                 proxy_BCA           = proxy_BCA,
                 proxy_CATE          = proxy_CATE,
                 membership          = membership,
                 X1_control          = X1_control,
                 diff                = diff,
                 vcov_control        = vcov_control,
                 monotonize          = monotonize,
                 external_weights    = external_weights,
                 significance_level  = significance_level)

} # FUN


# helper function that skips the input checks
GATES_NoChecks <- function(D, Y,
                           propensity_scores,
                           proxy_BCA,
                           proxy_CATE,
                           membership,
                           HT                  = FALSE,
                           X1_control          = setup_X1(),
                           diff                = setup_diff(),
                           vcov_control        = setup_vcov(),
                           monotonize          = TRUE,
                           external_weights    = NULL,
                           significance_level  = 0.05){

  # fit model according to strategy 1 or 2 in the paper
  do.call(what = get(ifelse(HT, "GATES.HT", "GATES.classic")),
          args = list(D = D, Y = Y,
                      propensity_scores   = propensity_scores,
                      proxy_BCA           = proxy_BCA,
                      proxy_CATE          = proxy_CATE,
                      membership          = membership,
                      X1_control          = X1_control,
                      diff                = diff,
                      vcov_control        = vcov_control,
                      monotonize          = monotonize,
                      external_weights    = external_weights,
                      significance_level  = significance_level))

} # FUN


# helper function for case when there is no HT transformation used. Wrapped by function "GATES"
GATES.classic <- function(D, Y,
                          propensity_scores,
                          proxy_BCA, proxy_CATE,
                          membership,
                          X1_control         = setup_X1(),
                          diff               = setup_diff(),
                          vcov_control       = setup_vcov(),
                          monotonize         = TRUE,
                          external_weights   = NULL,
                          significance_level = 0.05){

  # make the group membership a binary matrix
  groups <- 1 * membership

  # number of groups
  K <- ncol(groups)

  # prepare weights
  weights <- 1 / (propensity_scores * (1 - propensity_scores))

  # if external weights are supplied, include them in the weighting
  if(!is.null(external_weights))
  {
    weights <- weights * external_weights
  } # IF

  # prepare matrix X1
  X1     <- get.df.from.X1_control(functions.of.Z_mat = cbind(S = proxy_CATE,
                                                                B = proxy_BCA,
                                                                p = propensity_scores),
                                     X1_control = X1_control)

  # prepare covariate matrix
  X <- data.frame(X1,
                  (D - propensity_scores) * groups)
  colnames(X) <- c(colnames(X1), paste0("gamma.", 1:K))

  # fit weighted linear regression by OLS
  gates.obj <- stats::lm(Y ~., data = data.frame(Y, X), weights = weights)

  # get estimate of covariance matrix of the error terms
  vcov. <- get.vcov(x              = gates.obj,
                    vcov_control   = vcov_control)

  # extract the relevant coefficients
  coefficients                 <- lmtest::coeftest(gates.obj, vcov. = vcov.)

  # return
  return(structure(
    list(generic_targets = generic_targets_GATES(coeftest.object = coefficients,
                                                 K = K,
                                                 vcov = vcov.,
                                                 monotonize = monotonize,
                                                 significance_level = significance_level,
                                                 diff = diff),
         coefficients = coefficients,
         lm = gates.obj),
    class = "GATES"))

} # END FUN


# helper function for case when there is a HT transformation used. Wrapped by function "GATES"
GATES.HT <- function(D, Y,
                     propensity_scores,
                     proxy_BCA, proxy_CATE,
                     membership,
                     X1_control         = setup_X1(),
                     diff               = setup_diff(),
                     vcov_control       = setup_vcov(),
                     monotonize         = TRUE,
                     external_weights   = NULL,
                     significance_level = 0.05){

  # make the group membership a binary matrix
  groups <- 1 * membership

  # number of groups
  K <- ncol(groups)

  # HT transformation
  H <- (D - propensity_scores) / (propensity_scores * (1 - propensity_scores))

  # prepare matrix X1
  X1 <- get.df.from.X1_control(functions.of.Z_mat = cbind(S = proxy_CATE,
                                                            B = proxy_BCA,
                                                            p = propensity_scores),
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

  # prepare the external weights (if applicable)
  if(is.null(external_weights))
  {
    weights <- NULL
  } else{
    weights <- external_weights
  }

  # fit linear regression by OLS (no intercept!)
  gates.obj <- stats::lm(formula =  stats::as.formula(paste0("YH ~ ", paste0(colnames(X), collapse = " + "), " + 0")),
                data = data.frame(YH = Y*H, X), weights = weights)

  # get estimate of covariance matrix of the error terms
  vcov. <- get.vcov(x              = gates.obj,
                    vcov_control   = vcov_control)

  # extract the relevant coefficients
  coefficients                 <- lmtest::coeftest(gates.obj, vcov. = vcov.)

  # return
  return(structure(
    list(generic_targets = generic_targets_GATES(coeftest.object = coefficients,
                                                 K = K,
                                                 vcov = vcov.,
                                                 monotonize = monotonize,
                                                 significance_level = significance_level,
                                                 diff = diff),
         coefficients = coefficients,
         lm = gates.obj),
    class = "GATES"))

} # END FUN


# helper function to calculate the generic targets of BLP (vcov is estimate, not list)
generic_targets_GATES <- function(coeftest.object, K, vcov, monotonize,
                                  significance_level = 0.05,
                                  diff = diff){

  # only gamma covariances are relevant here
  gammanam <- paste0("gamma.", seq_len(K))
  vcov <- vcov[gammanam, gammanam]

  # extract controls
  subtract_from  <- diff$subtract_from
  subtracted <- diff$subtracted

  # extract coefficients
  coefficients.temp <- coeftest.object[paste0("gamma.", 1:K), 1:3]
  colnames(coefficients.temp) <- c("Estimate", "Std. Error", "z value")

  # get quantile
  z          <- stats::qnorm(1-significance_level/2)

  # compute relevant statistics
  p.right <- stats::pnorm(coefficients.temp[,"z value"], lower.tail = FALSE) # right p value: Pr(Z>z)
  p.left  <- stats::pnorm(coefficients.temp[,"z value"], lower.tail = TRUE)  # left p value: Pr(Z<z)
  ci.lo   <- coefficients.temp[,"Estimate"] - z * coefficients.temp[,"Std. Error"]
  ci.up   <- coefficients.temp[,"Estimate"] + z * coefficients.temp[,"Std. Error"]


  # monotonize GATES point estimates and confidence bounds if requested
  if(monotonize)
  {
    ## rearrange point estimates to be monotonic following Chernozhukov et al. (2009, https://doi.org/10.1093/biomet/asp030)
    order_gamma <- order(coefficients.temp[,"Estimate"], decreasing = FALSE)
    nam <- rownames(coefficients.temp)
    coefficients.temp <- coefficients.temp[order_gamma,]
    rownames(coefficients.temp) <- nam
    p.left <-  structure(p.left[order_gamma], names = nam)
    p.right <- structure(p.right[order_gamma], names = nam)

    ## same for confidence bounds; see paragraph 3 Chernozhukov et al. (2009, https://doi.org/10.1093/biomet/asp030)
    ci.lo <- structure(sort(ci.lo, decreasing = FALSE), names = nam)
    ci.up <- structure(sort(ci.up, decreasing = FALSE), names = nam)

    ## adjust order in covariance matrix accordingly
    vcov_gamma <- vcov[order_gamma,order_gamma]
    colnames(vcov_gamma) <- rownames(vcov_gamma) <- gammanam
    vcov <- vcov_gamma
  }


  # prepare generic targets of the gammas
  generic_targets <- cbind(coefficients.temp, ci.lo, ci.up, p.left, p.right)
  colnames(generic_targets) <- c("Estimate", "Std. Error", "z value",
                                 "CB lower", "CB upper", "Pr(<z)", "Pr(>z)")
  generic_targets <- generic_targets[,c("Estimate", "CB lower", "CB upper",
                                        "Std. Error", "z value", "Pr(<z)", "Pr(>z)")]

  # GATES differences: estimate mean and variance of a difference
  if(subtract_from == "least"){

    diff. <- generic_targets["gamma.1", "Estimate"] -
      generic_targets[paste0("gamma.", subtracted), "Estimate"]

    diff.se <- as.numeric(
      sqrt(vcov["gamma.1", "gamma.1"] +
                      diag(vcov[paste0("gamma.", subtracted),
                                paste0("gamma.", subtracted), drop = FALSE]) -
                      2 * vcov["gamma.1", paste0("gamma.", subtracted)]))

    nam <- paste0("gamma.1-gamma.", subtracted)

  } else{

    diff. <- generic_targets[paste0("gamma.", K), "Estimate"] -
      generic_targets[paste0("gamma.", subtracted), "Estimate"]

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
  p.right <- stats::pnorm(zstat, lower.tail = FALSE) # right p value: Pr(Z>z)
  p.left  <- stats::pnorm(zstat, lower.tail = TRUE)  # left p value: Pr(Z<z)

  diff.mat <- cbind(diff., ci.lo, ci.up, diff.se, zstat, p.left, p.right)
  rownames(diff.mat) <- nam

  # return final matrix
  rbind(generic_targets, diff.mat)

} # FUN
