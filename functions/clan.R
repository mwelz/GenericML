#' Estimates the CLAN parameters in the main sample
#' 
#' @param Z.clan.main.sample a matrix with _|M|_ rows. Each column represents a variable for which CLAN shall be performed.
#' @param group.membership.main.sample a logical matrix with _M_ rows that indicate 
#' the group memberships (such a matrix is returned by the function quantile.group())
#' @param equal.group.variances logical. If TRUE, the the two within-group variances of the most and least affected group are assumed to be equal. Default is FALSE.
#' @param differences.control a list with two elements called 'group.to.subtract.from' and 'groups.to.be.subtracted'. The first element ('group.to.subtract.from') denotes what shall be the base group to subtract from in CLAN; either "most" or "least". The second element ('groups.to.be.subtracted') are the groups to be subtracted from 'group.to.subtract.from', which is a subset of {1,...,K}, where K equals the number of groups.
#' @return The two CLAN parameters ("most" affected and "least" affected) for each variable in Z.clan.main.sample
#' 
#' @export
CLAN <- function(Z.clan.main.sample, 
                 group.membership.main.sample, 
                 equal.group.variances = FALSE,
                 differences.control = list(group.to.subtract.from = "most",
                                            groups.to.be.subtracted = 1),
                 significance.level = 0.05){
  
  # extract controls
  group.to.subtract.from  <- differences.control$group.to.subtract.from
  groups.to.be.subtracted <- differences.control$groups.to.be.subtracted

  K <- ncol(group.membership.main.sample)
  group.base <- ifelse(group.to.subtract.from == "least", 1, K)
  
  # initialize
  generic.targets   <- rep(list(NULL), ncol(Z.clan.main.sample))
  clan.coefficients <- matrix(NA_real_, 
                              nrow =  K + length(groups.to.be.subtracted), 
                              ncol = ncol(Z.clan.main.sample))
  z                 <- qnorm(1-significance.level/2) # the quantile
  
  # loop over the CLAN variables
  for(j in 1:ncol(Z.clan.main.sample)){
    
    # initialize matrix
    out.mat <- matrix(NA_real_,
                      nrow = K + length(groups.to.be.subtracted), 
                      ncol = 7) 
    ct <- 1 # initialize counter
    
    
    ### 1. get summary statistics for most and least affected group ----
    for(k in 1:K){
      
      ttest.deltak <- stats::t.test(Z.clan.main.sample[group.membership.main.sample[, k], j])
      ci.lo        <- ttest.deltak$estimate - z * ttest.deltak$stderr 
      ci.up        <- ttest.deltak$estimate + z * ttest.deltak$stderr 
      p.right      <- pnorm(ttest.deltak$statistic, lower.tail = FALSE) # right p-value: Pr(Z>z)
      p.left       <- pnorm(ttest.deltak$statistic, lower.tail = TRUE)  # left p-value: Pr(Z<z)
      out.mat[ct,] <- c(ttest.deltak$estimate, ci.lo, ci.up,
                        ttest.deltak$stderr, ttest.deltak$statistic, p.left, p.right)
      ct           <- ct + 1 # update counter

    } # FOR
    
    
    
    ### 2. get summary statistics for differences ----
    
    for(k in groups.to.be.subtracted){
      ttest.diff   <- stats::t.test(x = Z.clan.main.sample[group.membership.main.sample[, group.base], j], 
                                    y = Z.clan.main.sample[group.membership.main.sample[, k], j], 
                                    var.equal = equal.group.variances) # 2-sample t-test
      
      diff         <- ifelse(group.base == 1, 
                             out.mat[1,1] - out.mat[k,1], 
                             out.mat[K,1] - out.mat[k,1])
      diff.se      <- ttest.diff$stderr
      ci.lo        <- diff - z * diff.se
      ci.up        <- diff + z * diff.se
      z.diff       <- ttest.diff$statistic
      p.right      <- pnorm(z.diff, lower.tail = FALSE) # right p-value: Pr(Z>z)
      p.left       <- pnorm(z.diff, lower.tail = TRUE)  # left p-value: Pr(Z<z)
      out.mat[ct,] <- c(diff, ci.lo, ci.up, diff.se, z.diff, p.left, p.right)
      ct           <- ct + 1 # update counter
      
    } # FOR

    colnames(out.mat)     <- c("Estimate", "CB lower", "CB upper", 
                               "Std. Error", "z value", "Pr(<z)", "Pr(>z)")
    rownames(out.mat)     <- c(paste0("delta.", 1:K),
                               paste0(
                               "delta.", group.base, "-",
                               "delta.", groups.to.be.subtracted))
    # "delta.K-delta.1")
    generic.targets[[j]]  <- out.mat
    clan.coefficients[,j] <- out.mat[,1] 
    
  } # END FOR
  
  names(generic.targets) <- colnames(clan.coefficients) <- colnames(Z.clan.main.sample)
  rownames(clan.coefficients) <- rownames(out.mat)
  
  return(list(clan.coefficients = clan.coefficients,
              generic.targets   = generic.targets))
  
} # END FUN
