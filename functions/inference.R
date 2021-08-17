#' Estimates the CLAN parameters in the main sample
#' 
#' @param Z.clan.main.sample a matrix with _|M|_ rows. Each column represents a variable for which CLAN shall be performed.
#' @param group.membership.main.sample a logical matrix with _M_ rows that indicate 
#' the group memberships (such a matrix is returned by the function quantile.group())
#' @param equal.group.variances logical. If TRUE, the the two within-group variances of the most and least affected group are assumed to be equal. Default is FALSE.
#' @return The two CLAN parameters ("most" affected and "least" affected) for each variable in Z.clan.main.sample
#' 
#' @export
CLAN <- function(Z.clan.main.sample, 
                 group.membership.main.sample, 
                 equal.group.variances = FALSE,
                 significance.level = 0.05){
  
  K <- ncol(group.membership.main.sample)
  
  # initialize
  generic.targets   <- list()
  clan.coefficients <- matrix(NA_real_, nrow = 3, ncol = ncol(Z.clan.main.sample))
  
  # loop over the CLAN variables
  for(j in 1:ncol(Z.clan.main.sample)){
    
    # initialize matrix
    out.mat <- matrix(NA_real_, nrow = 3, ncol = 7)
    
    # get summary statistics for least affected group
    ttest.delta1 <- stats::t.test(Z.clan.main.sample[group.membership.main.sample[, 1], j])
    ci.lo        <- ttest.delta1$estimate - qnorm(1-significance.level/2) * ttest.delta1$stderr 
    ci.up        <- ttest.delta1$estimate + qnorm(1-significance.level/2) * ttest.delta1$stderr 
    p.right      <- pnorm(ttest.delta1$statistic, lower.tail = FALSE) # right p-value: Pr(Z>z)
    p.left       <- pnorm(ttest.delta1$statistic, lower.tail = TRUE)  # left p-value: Pr(Z<z)
    out.mat[1,]  <- c(ttest.delta1$estimate, ci.lo, ci.up,
                      ttest.delta1$stderr, ttest.delta1$statistic, p.left, p.right)
    
    # get summary statistics for most affected group
    ttest.deltaK <- stats::t.test(Z.clan.main.sample[group.membership.main.sample[, K], j])
    ci.lo        <- ttest.deltaK$estimate - qnorm(1-significance.level/2) * ttest.deltaK$stderr 
    ci.up        <- ttest.deltaK$estimate + qnorm(1-significance.level/2) * ttest.deltaK$stderr 
    p.right      <- pnorm(ttest.deltaK$statistic, lower.tail = FALSE) # right p-value: Pr(Z>z)
    p.left       <- pnorm(ttest.deltaK$statistic, lower.tail = TRUE)  # left p-value: Pr(Z<z)
    out.mat[2,]  <- c(ttest.deltaK$estimate, ci.lo, ci.up,
                      ttest.deltaK$stderr, ttest.deltaK$statistic, p.left, p.right)
    
    # get summary statistics for difference between most and least affected group
    ttest.diff <- stats::t.test(x = Z.clan.main.sample[group.membership.main.sample[, K], j], 
                                y = Z.clan.main.sample[group.membership.main.sample[, 1], j], 
                                var.equal = equal.group.variances) # 2-sample t-test
    
    diff        <- ttest.deltaK$estimate - ttest.delta1$estimate # difference in means
    diff.se     <- ttest.diff$stderr
    ci.lo       <- diff - qnorm(1-significance.level/2) * diff.se
    ci.up       <- diff + qnorm(1-significance.level/2) * diff.se
    z.diff      <- ttest.diff$statistic
    p.right     <- pnorm(z.diff, lower.tail = FALSE) # right p-value: Pr(Z>z)
    p.left      <- pnorm(z.diff, lower.tail = TRUE)  # left p-value: Pr(Z<z)
    out.mat[3,] <- c(diff, ci.lo, ci.up, diff.se, z.diff, p.left, p.right)
    
    colnames(out.mat)     <- c("Estimate", "CB lower", "CB upper", 
                               "Std. Error", "z value", "Pr(<z)", "Pr(>z)")
    rownames(out.mat)     <- c("delta.1", "delta.K", "delta.K-delta.1")
    generic.targets[[j]]  <- out.mat
    clan.coefficients[,j] <- out.mat[,1] 
    
  } # END FOR
  
  names(generic.targets) <- colnames(clan.coefficients) <- colnames(Z.clan.main.sample)
  rownames(clan.coefficients) <- c("delta.1", "delta.K", "delta.K-delta.1")
  
  return(list(clan.coefficients = clan.coefficients,
              generic.targets   = generic.targets))
  
} # END FUN



#' returns the two parameters that are used to find the best ML method
#' 
#' @param BLP.obj an object as returned by CLAN()
#' @param GATES.obj an object as returned by get.BLP.parameters()
#' @param proxy.cate.main.sample Proxy CATE estimators for the main sample
#' @param group.membership.main.sample a logical matrix with _M_ rows that indicate 
#' the group memberships (such a matrix is returned by the function quantile.group())
#' @return lambda and lambda.bar parameters
#' 
#' @export
best.ml.method.parameters <- function(BLP.obj,
                                      GATES.obj, 
                                      proxy.cate.main.sample, 
                                      group.membership.main.sample){
  
  return(list(lambda = as.numeric(BLP.obj$blp.coefficients["beta.2"]^2 * var(proxy.cate.main.sample)),
              lambda.bar = as.numeric(colMeans(group.membership.main.sample) %*%  GATES.obj$gates.coefficients^2)))
  
} # END FUN