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
#' @param X1.variables_BLP a list controlling the variables that shall be used in the matrix X1 for the BLP regression. The first element of the list, functions_of_Z, needs to be a subset of c("S", "B", "p"), where "p" corresponds to the propensity scores (default is "B"). The seconds element, custom_covariates, is an optional matrix/data frame of custom covariates that shall be included in X1 (default is NULL). The third element, fixed_effects, is a vector of integers, strings, or a factor thereof that indicates group membership of the observations: For each group, a fixed effect will be added (default is NULL). Note that in the final matrix X1, a constant 1 will be silently included so that the regression model has an intercept.
#' @param X1.variables_GATES a list controlling the variables that shall be used in the matrix X1 for the GATES regression. The first element of the list, functions_of_Z, needs to be a subset of c("S", "B", "p"), where "p" corresponds to the propensity scores (default is "B"). The seconds element, custom_covariates, is an optional matrix/data frame of custom covariates that shall be included in X1 (default is NULL). The third element, fixed_effects, is a vector of integers, strings, or a factor thereof that indicates group membership of the observations: For each group, a fixed effect will be added (default is NULL). Note that in the final matrix X1, a constant 1 will be silently included if no HT transformation is applied so that the regression model has an intercept.
#' @param quantile.cutoffs The cutoff points of quantiles that shall be used for GATES grouping. Default is `c(0.25, 0.5, 0.75)`, which corresponds to the quartiles.
#' @param vcov.control_BLP a list with two elements called 'estimator' and 'arguments'. The argument 'estimator' is a string specifying the covariance matrix estimator to be used in the BLP regression; specifies a covariance estimator function in the sandwich package (https://cran.r-project.org/web/packages/sandwich/sandwich.pdf). Recommended estimators are "vcovBS", "vcovCL", "vcovHAC", and "vcovHC". Default is 'vcovHC'. The element 'arguments' is a list of arguments that shall be passed to the function specified in the element 'estimator'. Default leads to the (homoskedastic) ordinary least squares covariance matrix estimator. See the reference manual of the sandwich package for details (https://cran.r-project.org/web/packages/sandwich/vignettes/sandwich.pdf).
#' @param vcov.control_GATES same as vcov.control_BLP, just for GATES regression
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
                      X1.variables_BLP           = list(functions_of_Z = c("B"),
                                                        custom_covariates = NULL,
                                                        fixed_effects = NULL),
                      X1.variables_GATES         = list(functions_of_Z = c("B"),
                                                        custom_covariates = NULL,
                                                        fixed_effects = NULL),
                      quantile.cutoffs = c(0.25, 0.5, 0.75),
                      vcov.control_BLP           = list(estimator = "vcovHC",
                                                        arguments = list(type = "const")),
                      vcov.control_GATES         = list(estimator = "vcovHC",
                                                        arguments = list(type = "const")),
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
                               X1.variables_BLP           = X1.variables_BLP,
                               X1.variables_GATES         = X1.variables_GATES,
                               HT.transformation          = HT.transformation,
                               vcov.control_BLP           = vcov.control_BLP,
                               vcov.control_GATES         = vcov.control_GATES,
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
              propensity.scores = list(estimates = propensity.scores,
                                       mlr3.objects = propensity.scores.obj$mlr3.objects),
              genericML.by.split = gen.ml.different.learners$genericML.by.split,
              splits = gen.ml.different.learners$splits, 
              arguments = list(quantile.cutoffs = c(0.25, 0.5, 0.75),
                               proportion.in.main.set = 0.5, 
                               significance.level = 0.05,
                               learners.genericML = learners.genericML,
                               learner.propensity.score = learner.propensity.score)))
  
} # FUN
