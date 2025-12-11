#' Calculate lower and upper median
#'
#' Calculates the lower and and median of a vector as proposed in Comment 4.2 in the paper.
#'
#' @param x A numeric vector.
#'
#' @return
#' A list with the upper, lower, and usual median (where the latter is the average of the former two).
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

  # get lower median and upper median
  lower <- med_lo(x)
  upper <- med_up(x)

  return(list(lower_median = lower,
              upper_median = upper,
              median = mean(c(lower, upper))))

} # END FUN


## lower median
med_lo <- function(x) stats::quantile(x, probs = 0.5, type = 1, names = FALSE)

## upper median
med_up <- function(x)
{
  x_rev <- -x # reverse order
  -stats::quantile(x_rev, probs = 0.5, type = 1, names = FALSE) # undo reversing
}



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

  ## number of groups
  num_groups <- length(cutoffs) + 1L

  # to have non-zero within-group variation, we require at least 2 observations per group
  n              <- length(x)
  group_size_min <- 2L

  ### obtain the grouping
  # we require a while loop here because grouping on the raw x might be illegal, that is, the groups' sizes do not correspond to what one would expect under a continuous variable (see above). Such violations can happen when sufficiently large subsets of x have zero variation. To overcome this problem (if it is present), we add tiny random noise to x and repeat the grouping until a legal grouping is found.
  # NB: this problem was first spotted by Lukas Kitzmueller, who proposed the original bugfix (which we have adapted since). Many thanks to Lukas!
  grouping_unifinished <- TRUE
  ct <- 0L

  while(grouping_unifinished)
  {
    # get empirical quantiles
    q <- stats::quantile(x, cutoffs)

    # get names of breaks
    qnam <- breaks_format(breaks = q, dig.lab = 3L)

    # initialize out matrix and helper objects
    out <- matrix(NA, n, num_groups)
    legal_grouping <- rep(TRUE, num_groups)
    groupnam <- rep(NA_character_, num_groups)

    for(k in seq_len(num_groups))
    {
      ## get the grouping
      if(k == 1L)
      {
        bool_k <- x < q[k]
        groupnam[k] <- paste0("(-Inf, ", qnam[k], ")")

      } else if(k == num_groups)
      {
        bool_k <- x >= q[k-1L]
        groupnam[k] <- paste0("[", qnam[k-1L], ", Inf)")
      } else
      {
        bool_k <- q[k-1L] <= x & x < q[k]
        groupnam[k] <- paste0("[", qnam[k-1L], ", ",  qnam[k], ")")
      } # IF

      ## check if grouping is legal
      # a grouping is illegal if group doesn't have minimum size
      # this can happen if there is very little variation in x
      size_k <- sum(bool_k)
      if(size_k < group_size_min)
      {
        legal_grouping[k] <- FALSE
      } # IF

      ## store the candidate grouping
      out[,k] <- bool_k
    } # FOR k

    ## column naming
    colnames(out) <- groupnam

    ## check if the grouping is complete and legal
    if(all(legal_grouping))
    {
      grouping_unifinished <- FALSE # legal grouping -> will break the while loop
    } else
    {
      ## illegal grouping: induce tiny noise to increase variation
      x <- x + stats::rnorm(n, mean = 0.0, sd = 0.001)
    } # IF

    ## update counter
    ct <- ct + 1L

    if(ct > 3L)
    {
      stop("The specified quantile cutoffs do not allow for a grouping that results in groups with nonzero within-group variation. Increase the expected group size through the quantile cutoffs argument.")
    } # IF
  } # WHILE

  # return
  return(structure(out, type = "quantile_group"))

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
