#' Performs generic ML 
#' 
#' @param Z A matrix or data frame of the covariates.
#' @param D A binary vector of treatment assignment. 
#' @param Y The response vector.
#' @param learner.propensity.score The machine learner to be used for estimating the propensity scores. Either `'elastic.net'`, `'random.forest'`, or `'tree'`. Can alternatively be specified by using `mlr3` syntax, for example `'mlr3::lrn("ranger", num.trees = 500)'`. See https://mlr3learners.mlr-org.com for a list of `mlr3` learners. Default is `'elastic.net'`.
#' @param learners.genericML A vector of strings specifying the machine learners to be used for estimating the BCA and CATE. Either `'elastic.net'`, `'random.forest'`, or `'tree'`. Can alternatively be specified by using `mlr3` syntax, for example `'mlr3::lrn("ranger", num.trees = 500)'`. See https://mlr3learners.mlr-org.com for a list of `mlr3` learners.
#' @param num.splits number of sample splits. Default is 100.
#' @param Z.clan A matrix of variables that shall be considered for the CLAN. If `NULL` (default), then `Z.clan = Z`, i.e. CLAN is performed for all variables in `Z`.
#' @param quantile.cutoffs The cutoff points of quantiles that shall be used for GATES grouping. Default is `c(0.25, 0.5, 0.75)`, which corresponds to the quartiles.
#' @param vcov.type_BLP a character string specifying the estimation type of the error covariance matrix in BLP. See sandwich::vcovHC for details. Default is "const" (for homoskedasticity)
#' @param vcov.type_GATES a character string specifying the estimation type of the error covariance matrix in GATES. See sandwich::vcovHC for details. Default is "const" (for homoskedasticity)
#' @param equal.group.variances_CLAN logical. If TRUE, the the two within-group variances of the most and least affected group in CLAN are assumed to be equal. Default is FALSE.
#' @param proportion.in.main.set proportion of samples that shall be in main set. Default is 0.5.
#' @param significance.level significance level for VEIN. Default is 0.05.
#' @param store.learners Logical. If TRUE, all intermediate results of the learners will be stored. Default is FALSE.
#' @param store.splits Logical. If `TRUE`, information on the sample splits will be stored. Default is `FALSE`.
#' 
#' @export
GenericML <- function(Z, D, Y, 
                      learner.propensity.score = "elastic.net", 
                      learners.genericML,
                      num.splits = 100,
                      Z.clan = NULL,
                      quantile.cutoffs = c(0.25, 0.5, 0.75),
                      vcov.type_BLP = "const",
                      vcov.type_GATES = "const",
                      equal.group.variances_CLAN = FALSE,
                      proportion.in.main.set = 0.5, 
                      significance.level = 0.05,
                      store.learners = FALSE,
                      store.splits = FALSE){
  
  ### step 1: compute propensity scores ----
  propensity.scores.obj <- propensity.score(Z = Z, D = D, 
                                            learner = make.mlr3.string(learner.propensity.score, 
                                                                       regr = FALSE))
  propensity.scores     <- propensity.scores.obj$propensity.scores
  
  ### step 2: for each ML method, do the generic ML analysis ----
  
  gen.ml.different.learners <- 
    generic.ml.across.learners(Z = Z, D = D, Y = Y, 
                               propensity.scores = propensity.scores, 
                               learners = learners.genericML, 
                               num.splits = num.splits,
                               Z.clan = Z.clan, 
                               vcov.type_BLP = vcov.type_BLP,
                               vcov.type_GATES = vcov.type_GATES,
                               equal.group.variances_CLAN = equal.group.variances_CLAN,
                               proportion.in.main.set = proportion.in.main.set, 
                               quantile.cutoffs = quantile.cutoffs,
                               significance.level = significance.level,
                               store.learners = store.learners,
                               store.splits = store.splits)
  
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
