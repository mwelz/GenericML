
# get lower and upper median as in comment 4.2 in Chernozhukov et al. (2021)
Med <- function(x){
  
  # get the empirical CDF of X
  ecdf.x <- ecdf(x)
  
  # evaluate the ensuing probabilities
  F.x <- ecdf.x(x)
  
  # get lower median and upper median
  lower       <- min(x[F.x >= 0.5])
  upper.array <- x[(1 - F.x) >= 0.5]
  upper       <- ifelse(length(upper.array) == 0, 
                        lower,
                        max(upper.array)) # account for case where upper.array is empty
  
  return(list(lower.median = lower, 
              upper.median = upper, 
              Med = mean(c(lower, upper))))
  
} # END FUN








#' Partition a vector x into groups based on its quantiles
#' 
#' @param x the vector to be partitioned
#' @param cutoffs the quantile cutoffs for the partition. Default are the quartiles.
#' @param quantile.nam logical. Shall the cutoff values be included?
quantile.group <- function(x,
                           cutoffs = c(0.25, 0.5, 0.75),
                           quantile.nam = TRUE){
  # cutoffs are the quantile cutoffs (like c(0.25, 0.5, 0.75))
  q         <- quantile(x, cutoffs)
  q         <- c(-Inf, q, Inf)
  groups    <- as.character(cut(x, breaks = q, include.lowest = TRUE, right = FALSE, dig.lab = 3))
  group.nam <- unique(groups)
  group.nam <- group.nam[order(
    as.numeric(substr(sub("\\,.*", "", group.nam), 2, stop = 1e8L)), 
    decreasing = FALSE)] # ensure the order is correct
  group.mat <- matrix(NA, length(x), length(group.nam))
  nam       <- rep(NA, length(group.nam))
  
  for(j in 1:length(group.nam)){
    if(j == 1){
      nam[j] <- paste0("<", 100*cutoffs[j], "% quantile")
    } else if (j == length(group.nam)){
      nam[j] <- paste0(">=", 100*cutoffs[j-1], "% quantile")
    } else{
      nam[j] <- paste0("[", 100*cutoffs[j-1], ",", 100*cutoffs[j], ")% quantile")
    }
    group.mat[,j] <- groups == group.nam[j]
  }
  
  if(quantile.nam){
    colnames(group.mat) <- nam
  } else{
    colnames(group.mat) <- group.nam
  }
  return(group.mat)
} # FUN



#' Performs generic ML for a given learning technique (with only one split of the data)
#' 
#' @param Z a matrix or data frame of covariates
#' @param D a binary vector of treatment status of length 
#' @param Y a vector of responses of length
#' @param propensity.scores a vector of propensity scores
#' @param learner The machine learner that shall be used
#' @param M.set main set
#' @param A.set auxiliary set
#' @param vcov.type_BLP a character string specifying the estimation type of the error covariance matrix in BLP. See sandwich::vcovHC for details. Default is "const" (for homoskedasticity)
#' @param vcov.type_GATES a character string specifying the estimation type of the error covariance matrix in GATES. See sandwich::vcovHC for details. Default is "const" (for homoskedasticity)
#' @param equal.group.variances_CLAN logical. If TRUE, the the two within-group variances of the most and least affected group in CLAN are assumed to be equal. Default is FALSE.

#' @param quantile.cutoffs Cutoff points of quantiles that shall be used for GATES grouping
#' 
#' TODO: instructions on how mlr3 input is supposed to work (needs to be a string!)
#' TODO: comments on CLAN: If there are categorical variables, apply one-hot-encoding to Z.clan. The interpretation then becomes: Is there a factor that is overproportionally present in the least or most affected group?
#' 
#' @export
get.generic.ml.for.given.learner <- function(Z, D, Y, 
                                             propensity.scores,
                                             learner = 'mlr3::lrn("cv_glmnet", s = "lambda.min")',
                                             M.set, A.set,
                                             Z.clan = NULL, 
                                             vcov.type_BLP = "const",
                                             vcov.type_GATES = "const",
                                             equal.group.variances_CLAN = FALSE,
                                             proportion.in.main.set = 0.5, 
                                             quantile.cutoffs = c(0.25, 0.5, 0.75),
                                             significance.level = 0.05){
  
  ### step 1: input checks ---- 
  if(is.null(Z.clan)) Z.clan <- Z # if no input provided, set it equal to Z
  
  ### step 2a: learn proxy predictors by using the auxiliary set ----
  
  # get the proxy baseline estimator for the main sample
  proxy.baseline.obj <- baseline.proxy.estimator(Z = Z, D = D, Y = Y, 
                                                 auxiliary.sample = A.set, 
                                                 learner = make.mlr3.string(learner, regr = TRUE))
  proxy.baseline     <- proxy.baseline.obj$baseline.predictions.main.sample
  
  # get the proxy estimator of the CATE for the main sample
  proxy.cate.obj <- 
    CATE.proxy.estimator(Z = Z, D = D, Y = Y,
                         auxiliary.sample = A.set, 
                         learner = make.mlr3.string(learner, regr = TRUE),
                         proxy.baseline.estimates = proxy.baseline.obj$baseline.predictions.full.sample)
  proxy.cate <- proxy.cate.obj$CATE.predictions.main.sample
  
  
  ### step 2b: estimate BLP parameters by OLS (TODO: HT transformation!) ----
  blp.obj <- BLP.classic(D = D[M.set], Y = Y[M.set],
                         propensity.scores = propensity.scores[M.set],
                         proxy.baseline = proxy.baseline, 
                         proxy.cate = proxy.cate, 
                         vcov.type = vcov.type_BLP,
                         significance.level = significance.level)
  
  
  ### step 2c: estimate GATES parameters by OLS (TODO: HT transformation!) ----
  # group the proxy estimators for the CATE in the main sample by quantiles
  group.membership.main.sample <- quantile.group(proxy.cate, 
                                                 cutoffs = quantile.cutoffs, 
                                                 quantile.nam = TRUE) 
  
  gates.obj <- GATES.classic(D = D[M.set], Y = Y[M.set],
                             propensity.scores = propensity.scores[M.set],
                             group.membership.main.sample = group.membership.main.sample, 
                             proxy.baseline = proxy.baseline, proxy.cate = proxy.cate,
                             vcov.type = vcov.type_GATES,
                             significance.level = significance.level)
  
  
  ### step 2d: estimate CLAN parameters in the main sample
  clan.obj <- CLAN(Z.clan.main.sample = Z.clan[M.set,], 
                   group.membership.main.sample = group.membership.main.sample,
                   equal.group.variances = equal.group.variances_CLAN,
                   significance.level = significance.level)
  
  
  ### step 2e: get parameters over which we maximize to find the "best" ML method ----
  best.obj <- best.ml.method.parameters(BLP.obj = blp.obj, 
                                        GATES.obj = gates.obj, 
                                        proxy.cate.main.sample = proxy.cate,
                                        group.membership.main.sample = group.membership.main.sample)
  
  ### organize output in a list ----
  return(list(BLP = blp.obj, 
              GATES = gates.obj, 
              CLAN = clan.obj,
              best = best.obj,
              CATE.proxy = proxy.cate.obj,
              baseline.proxy = proxy.baseline,
              group.membership_M.set = group.membership.main.sample))
  
} # END FUN