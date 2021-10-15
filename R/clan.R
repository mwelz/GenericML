#' Estimates the CLAN parameters in the main sample
#'
#' @param Z_CLAN.main.sample a matrix with _|M|_ rows. Each column represents a variable for which CLAN shall be performed.
#' @param group.membership.main.sample a logical matrix with _M_ rows that indicate
#' the group memberships (such a matrix is returned by the function quantile.group())
#' @param equal.group.variances logical. If TRUE, the the two within-group variances of the most and least affected group are assumed to be equal. Default is FALSE.
#' @param diff Controls the generic targets of CLAN. See the documentation of \code{\link{setup_diff}}.
#' @param significance_level Significance level. Default is 0.05.
#' @return The two CLAN parameters ("most" affected and "least" affected) for each variable in Z_CLAN.main.sample
#'
#' @export
CLAN <- function(Z_CLAN.main.sample,
                 group.membership.main.sample,
                 equal.group.variances = FALSE,
                 diff = setup_diff(),
                 significance_level = 0.05){

  # input checks
  InputChecks_group.membership(group.membership.main.sample)
  InputChecks_equal.length2(Z_CLAN.main.sample, group.membership.main.sample)
  InputChecks_diff(diff, K = ncol(group.membership.main.sample))

  # run main function
  CLAN_NoChecks(Z_CLAN.main.sample = Z_CLAN.main.sample,
                group.membership.main.sample = group.membership.main.sample,
                equal.group.variances = equal.group.variances,
                diff = diff,
                significance_level = significance_level)

} # END FUN



# performs CLAN without calling the input checks
CLAN_NoChecks <- function(Z_CLAN.main.sample,
                          group.membership.main.sample,
                          equal.group.variances = FALSE,
                          diff = setup_diff(),
                          significance_level = 0.05){

  # extract controls
  subtract_from  <- diff$subtract_from
  subtracted <- diff$subtracted

  K <- ncol(group.membership.main.sample)
  group.base <- ifelse(subtract_from == "least", 1, K)

  # initialize
  generic.targets   <- rep(list(NULL), ncol(Z_CLAN.main.sample))
  clan.coefficients <- matrix(NA_real_,
                              nrow =  K + length(subtracted),
                              ncol = ncol(Z_CLAN.main.sample))
  z                 <- stats::qnorm(1-significance_level/2) # the quantile

  # loop over the CLAN variables
  for(j in 1:ncol(Z_CLAN.main.sample)){

    # initialize matrix
    out.mat <- matrix(NA_real_,
                      nrow = K + length(subtracted),
                      ncol = 7)
    ct <- 1 # initialize counter


    ### 1. get summary statistics for most and least affected group ----
    for(k in 1:K){

      if(stats::var(Z_CLAN.main.sample[group.membership.main.sample[, k], j]) == 0){

        # in case of zero variation, t.test() will throw an error. In this case, return uninformative out.mat[ct,]. NB: this bug has been spotted and fixed by Lucas Kitzmueller. All credits for this fix go to him!
        mean.estimate <- mean(Z_CLAN.main.sample[group.membership.main.sample[, k], j])
        out.mat[ct,]  <- c(mean.estimate, mean.estimate, mean.estimate,
                           0.0, 0.0, 0.5, 0.5)

      } else{

        ttest.deltak <- stats::t.test(Z_CLAN.main.sample[group.membership.main.sample[, k], j])
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

      x <- Z_CLAN.main.sample[group.membership.main.sample[, group.base], j]
      y <- Z_CLAN.main.sample[group.membership.main.sample[, k], j]

      if((stats::var(x) == stats::var(y)) & (stats::var(x) == 0)){

        # in case of zero variation, t.test() will throw an error. In this case, return uninformative out.mat[ct,]. NB: this bug has been spotted and fixed by Lucas Kitzmueller. All credits for this fix go to him!
        diff. <- ci.lo <- ci.up <- mean(x) - mean(y)
        diff.se <- z.diff      <- 0.0
        p.left <- p.right      <- 0.5

      } else{

        ttest.diff   <- stats::t.test(x = x, y = y,
                                      var.equal = equal.group.variances) # 2-sample t-test
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
    generic.targets[[j]]  <- out.mat
    clan.coefficients[,j] <- out.mat[,1]

  } # END FOR

  names(generic.targets) <- colnames(clan.coefficients) <- colnames(Z_CLAN.main.sample)
  rownames(clan.coefficients) <- rownames(out.mat)

  return(structure(
    list(clan.coefficients = clan.coefficients,
         generic.targets   = generic.targets), class = "CLAN"))

} # FUN
