#' Estimates the BLP parameters based on the main sample M. 
#' 
#' @param D a binary vector of treatment status of length _|M|_
#' @param Y a vector of responses of length _|M|_
#' @param proxy.baseline a vector of proxy baseline estimates of length _M_
#' @param proxy.cate a vector of proxy CATE estimates of length _M_
#' @param propensity.scores a vector of propensity scores of length _|M|_
#' @return BLP coefficients with inference statements
#' 
#' @export
#' 
#' TODO: implement same with HT transformation! 
BLP.classic <- function(D, Y, propensity.scores, 
                        proxy.baseline, proxy.cate, 
                        significance.level = 0.05){
  
  # prepare weights
  weights <- 1 / (propensity.scores * (1 - propensity.scores))
  
  # prepare covariate matrix
  X <- data.frame(B = proxy.baseline, 
                  beta.1 = D - propensity.scores, 
                  beta.2 = (D - propensity.scores) * (proxy.cate - mean(proxy.cate))) 
  
  # fit weighted linear regression by OLS
  blp.obj <- lm(Y ~., data = data.frame(Y, X), weights = weights)
  
  # extract coefficients
  coefficients     <- summary(blp.obj)$coefficients
  
  # generic targets
  coefficients.temp <- coefficients[c("beta.1", "beta.2"), 1:3]
  colnames(coefficients.temp) <- c("Estimate", "Std. Error", "z value")
  p.right <- pnorm(coefficients.temp[,"z value"], lower.tail = FALSE) # right p-value: Pr(Z>z)
  p.left  <- pnorm(coefficients.temp[,"z value"], lower.tail = TRUE)  # left p-value: Pr(Z<z)
  ci.lo   <- coefficients.temp[,"Estimate"] - qnorm(1-significance.level/2) * coefficients.temp[,"Std. Error"]
  ci.up   <- coefficients.temp[,"Estimate"] + qnorm(1-significance.level/2) * coefficients.temp[,"Std. Error"]
  generic.targets <- cbind(coefficients.temp, ci.lo, ci.up, p.left, p.right)
  colnames(generic.targets) <- c("Estimate", "Std. Error", "z value", 
                                 "CB lower", "CB upper", "Pr(<z)", "Pr(>z)")
  generic.targets <- generic.targets[,c("Estimate", "CB lower", "CB upper", 
                                        "Std. Error", "z value", "Pr(<z)", "Pr(>z)")]
  
  return(list(lm.obj = blp.obj, 
              blp.coefficients = blp.obj$coefficients[c("beta.1", "beta.2")],
              generic.targets = generic.targets,
              coefficients = coefficients))
  
} # END FUN



#' Estimates the GATES parameters based on the main sample M. 
#' 
#' @param D a binary vector of treatment status of length _|M|_
#' @param Y a vector of responses of length _|M|_
#' @param propensity.scores a vector of propensity scores of length _|M|_
#' @param proxy.baseline a vector of proxy baseline estimates of length _M_
#' @param proxy.cate a vector of proxy CATE estimates of length _M_
#' @param group.membership.main.sample a logical matrix with _M_ rows that indicate 
#' the group memberships (such a matrix is returned by the function quantile.group())
#' @return GATES coefficients 
#' 
#' @export
#' 
#' TODO: implement same with HT transformation! 
GATES.classic <- function(D, Y, 
                          propensity.scores, 
                          proxy.baseline, proxy.cate,
                          group.membership.main.sample,
                          significance.level = 0.05){
  
  # make the group membership a binary matrix
  groups <- 1 * group.membership.main.sample
  
  # prepare weights
  weights <- 1 / (propensity.scores * (1 - propensity.scores))
  
  # prepare covariate matrix
  X <- data.frame(B = proxy.baseline, 
                  (D - propensity.scores) * groups)
  colnames(X) <- c("B", paste0("gamma.", 1:ncol(groups)))
  
  # fit weighted linear regression by OLS
  gates.obj <- lm(Y ~., data = data.frame(Y, X), weights = weights)
  
  # extract coefficients
  coefficients                 <- summary(gates.obj)$coefficients
  gates.coefficients           <- coefficients[paste0("gamma.", 1:ncol(groups)), 1]
  gates.coefficients.quantiles <- colnames(groups)
  names(gates.coefficients.quantiles) <- paste0("gamma.", 1:ncol(groups))
  
  # generic targets
  coefficients.temp <- coefficients[-c(1,2), 1:3]
  colnames(coefficients.temp) <- c("Estimate", "Std. Error", "z value")
  p.right <- pnorm(coefficients.temp[,"z value"], lower.tail = FALSE) # right p-value: Pr(Z>z)
  p.left  <- pnorm(coefficients.temp[,"z value"], lower.tail = TRUE)  # left p-value: Pr(Z<z)
  ci.lo   <- coefficients.temp[,"Estimate"] - qnorm(1-significance.level/2) * coefficients.temp[,"Std. Error"]
  ci.up   <- coefficients.temp[,"Estimate"] + qnorm(1-significance.level/2) * coefficients.temp[,"Std. Error"]
  generic.targets <- cbind(coefficients.temp, ci.lo, ci.up, p.left, p.right)
  colnames(generic.targets) <- c("Estimate", "Std. Error", "z value", 
                                 "CB lower", "CB upper", "Pr(<z)", "Pr(>z)")
  generic.targets <- generic.targets[,c("Estimate", "CB lower", "CB upper", 
                                        "Std. Error", "z value", "Pr(<z)", "Pr(>z)")]
  
  # prepare generic target parameters for the difference
  covmat  <- stats::vcov(gates.obj)
  diff    <- coefficients[paste0("gamma.", ncol(groups)), "Estimate"] - 
    coefficients["gamma.1", "Estimate"]
  diff.se <- sqrt(covmat[paste0("gamma.", ncol(groups)), paste0("gamma.", ncol(groups))] +
                    covmat["gamma.1", "gamma.1"] - 2 * covmat[paste0("gamma.", ncol(groups)), "gamma.1"])
  ci.lo   <- diff - qnorm(1-significance.level/2) * diff.se
  ci.up   <- diff + qnorm(1-significance.level/2) * diff.se
  zstat   <- diff / diff.se
  p.right <- pnorm(zstat, lower.tail = FALSE) # right p-value: Pr(Z>z)
  p.left  <- pnorm(zstat, lower.tail = TRUE)  # left p-value: Pr(Z<z)
  
  return(list(lm.obj = gates.obj, 
              gates.coefficients = gates.coefficients,
              gates.coefficients.quantiles = gates.coefficients.quantiles,
              generic.targets = rbind(generic.targets, 
                                      matrix(c(diff, ci.lo, ci.up, diff.se, zstat, p.left, p.right), nrow = 1, 
                                             dimnames = list("gamma.K-gamma.1", NULL)) ),
              coefficients = coefficients))
  
} # END FUN



#' Estimates the CLAN parameters in the main sample
#' 
#' @param Z.clan.main.sample a matrix with _|M|_ rows. Each column represents a variable for which CLAN shall be performed.
#' @param group.membership.main.sample a logical matrix with _M_ rows that indicate 
#' the group memberships (such a matrix is returned by the function quantile.group())
#' @return The two CLAN parameters ("most" affected and "least" affected) for each variable in Z.clan.main.sample
#' 
#' @export
CLAN <- function(Z.clan.main.sample, 
                 group.membership.main.sample, 
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
                                var.equal = FALSE) # 2-sample Welch t-test
    
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