#' Performs generic ML 
#' 
#' @param Z A matrix or data frame of the covariates.
#' @param D A binary vector of treatment assignment. 
#' @param Y The response vector.
#' @param learner.propensity.score the estimator to be used. Either a numeric vector (which is then taken as estimates of the propensity scores) or a string specifying the estimator. The string must either be equal to 'constant' (estimates the propensity scores by mean(D)), 'elastic.net', 'random.forest', 'tree', or mlr3 syntax. Example for the latter: mlr3::lrn("classif.ranger", num.trees = 500) for a classification forest.
#' @param learners.genericML A vector of strings specifying the machine learners to be used for estimating the BCA and CATE. Either `'elastic.net'`, `'random.forest'`, or `'tree'`. Can alternatively be specified by using `mlr3` syntax, for example `'mlr3::lrn("ranger", num.trees = 500)'`. See https://mlr3learners.mlr-org.com for a list of `mlr3` learners.
#' @param num.splits number of sample splits. Default is 100.
#' @param Z.clan A matrix of variables that shall be considered for the CLAN. If `NULL` (default), then `Z.clan = Z`, i.e. CLAN is performed for all variables in `Z`.
#' @param HT.transformation logical. If TRUE, a HT transformation is applied in BLP and GATES. Default is FALSE.
#' @param X1.variables a character string specifying the variables in the matrix X1. Needs to be a subset of c("S", "B", "p"), where "p" corresponds to the propensity scores. Unless a HT transformation is employed in GATES, a constant 1 is silently included in X1 as well.
#' @param quantile.cutoffs The cutoff points of quantiles that shall be used for GATES grouping. Default is `c(0.25, 0.5, 0.75)`, which corresponds to the quartiles.
#' @param vcov.type_BLP a character string specifying the estimation type of the error covariance matrix in BLP. See sandwich::vcovHC for details. Default is "const" (for homoskedasticity)
#' @param vcov.type_GATES a character string specifying the estimation type of the error covariance matrix in GATES. See sandwich::vcovHC for details. Default is "const" (for homoskedasticity)
#' @param equal.group.variances_CLAN logical. If TRUE, the the two within-group variances of the most and least affected group in CLAN are assumed to be equal. Default is FALSE.
#' @param proportion.in.main.set proportion of samples that shall be in main set. Default is 0.5.
#' @param significance.level significance level for VEIN. Default is 0.05.
#' @param minimum.variation minimum variation of the predictions before random noise with distribution N(0, var(Y)/20) is added. Default is 1e-05.
#' @param store.learners Logical. If TRUE, all intermediate results of the learners will be stored. Default is FALSE.
#' @param store.splits Logical. If `TRUE`, information on the sample splits will be stored. Default is `FALSE`.
#' 
#' @export
GenericML <- function(Z, D, Y, 
                      learner.propensity.score = "constant", 
                      learners.genericML,
                      num.splits = 100,
                      Z.clan = NULL,
                      HT.transformation = FALSE,
                      X1.variables = c("B"),
                      quantile.cutoffs = c(0.25, 0.5, 0.75),
                      vcov.type_BLP = "const",
                      vcov.type_GATES = "const",
                      equal.group.variances_CLAN = FALSE,
                      proportion.in.main.set = 0.5, 
                      significance.level = 0.05,
                      minimum.variation = 1e-05,
                      store.learners = FALSE,
                      store.splits = FALSE){
  
  ### step 1: compute propensity scores ----
  propensity.scores.obj <- propensity.score(Z = Z, D = D, 
                                            estimator = learner.propensity.score)
  propensity.scores     <- propensity.scores.obj$propensity.scores
  
  ### step 2: for each ML method, do the generic ML analysis ----
  
  gen.ml.different.learners <- 
    generic.ml.across.learners(Z = Z, D = D, Y = Y, 
                               propensity.scores          = propensity.scores, 
                               learners                   = learners.genericML, 
                               num.splits                 = num.splits,
                               Z.clan                     = Z.clan, 
                               X1.variables               = X1.variables,
                               HT.transformation          = HT.transformation,
                               vcov.type_BLP              = vcov.type_BLP,
                               vcov.type_GATES            = vcov.type_GATES,
                               equal.group.variances_CLAN = equal.group.variances_CLAN,
                               proportion.in.main.set     = proportion.in.main.set, 
                               quantile.cutoffs           = quantile.cutoffs,
                               significance.level         = significance.level,
                               minimum.variation          = minimum.variation,
                               store.learners             = store.learners,
                               store.splits               = store.splits)
  
  # extract the best learners
  best.learners <- get.best.learners(gen.ml.different.learners$generic.targets)
  
  ### step 3: perform VEIN analysis ---- 
  vein <- VEIN(gen.ml.different.learners$generic.targets, best.learners)
  
  return(list(VEIN = vein,
              best.learners = best.learners,
              propensity.scores = propensity.scores,
              genericML.by.split = gen.ml.different.learners$genericML.by.split,
              splits = gen.ml.different.learners$splits, 
              arguments = list(quantile.cutoffs = c(0.25, 0.5, 0.75),
                               proportion.in.main.set = 0.5, 
                               significance.level = 0.05,
                               learners.genericML = learners.genericML,
                               learner.propensity.score = learner.propensity.score)))
  
} # FUN
