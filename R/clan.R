#' Performs CLAN
#'
#' Performs Classification Analysis (CLAN) on all variables in a design matrix.
#'
#' @param Z_CLAN A numeric matrix holding variables on which classification analysis (CLAN) shall be performed. CLAN will be performed on each column of the matrix.
#' @param membership A logical matrix that indicates the group membership of each observation in \code{Z_CLAN}. Needs to be of type \code{\link{quantile_group}}. Typically, the grouping is based on CATE estimates, which are for instance returned by \code{proxy_CATE}.
#' @param equal_variances If \code{TRUE}, then all within-group variances of the CLAN groups are assumed to be equal. Default is \code{FALSE}. This specification is required for heteroskedasticity-robust variance estimation on the difference of two CLAN generic targets (i.e. variance of the difference of two means). If \code{TRUE} (corresponds to homoskedasticity assumption), the pooled variance is used. If \code{FALSE} (heteroskedasticity), the variance of Welch's t-test is used.
#' @param diff Specifies the generic targets of CLAN. Must be an instance of \code{\link{setup_diff}}. See the documentation of \code{\link{setup_diff}} for details.
#' @param significance_level Significance level. Default is 0.05.
#'
#' @return An object of the class \code{CLAN}, consisting of the following components:
#' \describe{
#'   \item{\code{generic_targets}}{A list of result matrices for each variable in \code{Z_CLAN}. Each matrix contains inferential results on the CLAN generic targets.}
#'   \item{\code{coefficients}}{A matrix of point estimates of each CLAN generic target parameter.}
#'   }
#'
#' @seealso
#' \code{\link{quantile_group}},
#' \code{\link{setup_diff}}
#'
#' @references
#' Chernozhukov V., Demirer M., Duflo E., Fernández-Val I. (2020). \dQuote{Generic Machine Learning Inference on Heterogenous Treatment Effects in Randomized Experiments.} \emph{arXiv preprint arXiv:1712.04802}. URL: \url{https://arxiv.org/abs/1712.04802}.
#'
#' @examples
#' ## generate data
#' set.seed(1)
#' n  <- 150                              # number of observations
#' p  <- 5                                # number of covariates
#' Z_CLAN <- matrix(runif(n*p), n, p)     # design matrix to perform CLAN on
#' membership <- quantile_group(rnorm(n)) # group membership
#'
#' ## perform CLAN
#' CLAN(Z_CLAN, membership)
#'
#' @export
CLAN <- function(Z_CLAN,
                 membership,
                 equal_variances    = FALSE,
                 diff               = setup_diff(),
                 significance_level = 0.05){

  # input checks
  InputChecks_Z_CLAN(Z_CLAN)
  InputChecks_equal.length2(Z_CLAN, membership)
  stopifnot(is.numeric(significance_level) & length(significance_level) == 1)
  stopifnot(0.0 < significance_level & significance_level < 0.5)
  stopifnot(is.logical(equal_variances))
  InputChecks_group.membership(membership)
  InputChecks_diff(diff, K = ncol(membership))

  # assign variable names if there are none
  if(is.null(colnames(Z_CLAN))) colnames(Z_CLAN) <- paste0("V", 1:ncol(Z_CLAN))

  # run main function
  CLAN_NoChecks(Z_CLAN = Z_CLAN,
                membership = membership,
                equal_variances = equal_variances,
                diff = diff,
                significance_level = significance_level)

} # END FUN



# performs CLAN without calling the input checks
CLAN_NoChecks <- function(Z_CLAN,
                          membership,
                          equal_variances = FALSE,
                          diff = setup_diff(),
                          significance_level = 0.05){

  # extract controls
  subtract_from  <- diff$subtract_from
  subtracted <- diff$subtracted

  K <- ncol(membership)
  group.base <- ifelse(subtract_from == "least", 1, K)

  # initialize
  generic_targets   <- rep(list(NULL), ncol(Z_CLAN))
  clan.coefficients <- matrix(NA_real_,
                              nrow =  K + length(subtracted),
                              ncol = ncol(Z_CLAN))
  z                 <- stats::qnorm(1-significance_level/2) # the quantile

  # loop over the CLAN variables
  for(j in 1:ncol(Z_CLAN)){

    # initialize matrix
    out.mat <- matrix(NA_real_,
                      nrow = K + length(subtracted),
                      ncol = 7)
    ct <- 1 # initialize counter


    ### 1. get summary statistics for most and least affected group ----
    for(k in 1:K){

      if(stats::var(Z_CLAN[membership[, k], j]) == 0){

        # in case of zero variation, t.test() will throw an error. In this case, return uninformative out.mat[ct,]. NB: this bug has been spotted and fixed by Lucas Kitzmueller. All credits for this fix go to him!
        mean.estimate <- mean(Z_CLAN[membership[, k], j])
        out.mat[ct,]  <- c(mean.estimate, mean.estimate, mean.estimate,
                           0.0, 0.0, 0.5, 0.5)

      } else{

        ttest.deltak <- stats::t.test(Z_CLAN[membership[, k], j])
        ci.lo        <- ttest.deltak$estimate - z * ttest.deltak$stderr
        ci.up        <- ttest.deltak$estimate + z * ttest.deltak$stderr
        p.right      <- stats::pnorm(ttest.deltak$statistic, lower.tail = FALSE) # right p-value: Pr(Z>z)
        p.left       <- stats::pnorm(ttest.deltak$statistic, lower.tail = TRUE)  # left p-value: Pr(Z<z)
        out.mat[ct,] <- c(ttest.deltak$estimate, ci.lo, ci.up,
                          ttest.deltak$stderr, ttest.deltak$statistic, p.left, p.right)

      } # IF

      ct           <- ct + 1 # update counter

    } # FOR



    ### 2. get summary statistics for differences ----

    for(k in subtracted){

      x <- Z_CLAN[membership[, group.base], j]
      y <- Z_CLAN[membership[, k], j]

      if((stats::var(x) == stats::var(y)) & (stats::var(x) == 0)){

        # in case of zero variation, t.test() will throw an error. In this case, return uninformative out.mat[ct,]. NB: this bug has been spotted and fixed by Lucas Kitzmueller. All credits for this fix go to him!
        diff. <- ci.lo <- ci.up <- mean(x) - mean(y)
        diff.se <- z.diff      <- 0.0
        p.left <- p.right      <- 0.5

      } else{

        ttest.diff   <- stats::t.test(x = x, y = y,
                                      var.equal = equal_variances) # 2-sample t-test
        diff.        <- ifelse(group.base == 1,
                               out.mat[1,1] - out.mat[k,1],
                               out.mat[K,1] - out.mat[k,1])
        diff.se      <- ttest.diff$stderr
        ci.lo        <- diff. - z * diff.se
        ci.up        <- diff. + z * diff.se
        z.diff       <- ttest.diff$statistic
        p.right      <- stats::pnorm(z.diff, lower.tail = FALSE) # right p-value: Pr(Z>z)
        p.left       <- stats::pnorm(z.diff, lower.tail = TRUE)  # left p-value: Pr(Z<z)

      } # IF

      out.mat[ct,] <- c(diff., ci.lo, ci.up, diff.se, z.diff, p.left, p.right)
      ct           <- ct + 1 # update counter

    } # FOR

    colnames(out.mat)     <- c("Estimate", "CB lower", "CB upper",
                               "Std. Error", "z value", "Pr(<z)", "Pr(>z)")
    rownames(out.mat)     <- c(paste0("delta.", 1:K),
                               paste0(
                                 "delta.", group.base, "-",
                                 "delta.", subtracted))
    # "delta.K-delta.1")
    generic_targets[[j]]  <- out.mat
    clan.coefficients[,j] <- out.mat[,1]

  } # END FOR

  names(generic_targets) <- colnames(clan.coefficients) <- colnames(Z_CLAN)
  rownames(clan.coefficients) <- rownames(out.mat)

  return(structure(
    list(generic_targets   = generic_targets,
         coefficients      = clan.coefficients),
    class = "CLAN"))

} # FUN
