#' Generic Machine Learning Inference
#'
#' Performs generic machine learning inference as in Chernozhukov, Demirer, Duflo and Fernández-Val (2020). Link to working paper: \url{https://arxiv.org/abs/1712.04802}.
#'
#' @param Z A matrix or data frame of the covariates.
#' @param D A binary vector of treatment assignment.
#' @param Y The response vector.
#' @param learner.propensity.score The estimator for the propensity scores. Either a numeric vector (which is then taken as estimates of the propensity scores) or a string specifying the estimator. The string must either be equal to \code{'constant'} (estimates the propensity scores by \code{mean(D)}), \code{'elastic.net'}, \code{'random.forest'}, \code{'tree'}, or \code{mlr3} syntax. Note that in case of \code{mlr3} syntax, do \emph{not} specify if the learner is a regression learner or classification learner; Example: \code{'mlr3::lrn("ranger", num.trees = 500)'} for a random forest learner. Note that this is a string and the absence of the \code{classif.} or \code{regr.} keyords.
#' @param learners.genericML A vector of strings specifying the machine learners to be used for estimating the BCA and CATE. Either \code{'elastic.net'}, \code{'random.forest'}, or \code{'tree'}. Can alternatively be specified by using \code{mlr3} syntax \emph{without} specification if the learner is a regression learner or classification learner. Example: \code{'mlr3::lrn("ranger", num.trees = 500)'} for a random forest learner. Note that this is a string and the absence of the \code{classif.} or \code{regr.} keyords. See \url{https://mlr3learners.mlr-org.com} for a list of \code{mlr3} learners.
#' @param num.splits Number of sample splits. Default is 100.
#' @param Z_CLAN A matrix of variables that shall be considered for the CLAN. If \code{NULL} (default), then \code{Z_CLAN = Z}, i.e. CLAN is performed for all variables in \code{Z}.
#' @param HT.transformation Logical. If TRUE, a HT transformation is applied in BLP and GATES. Default is FALSE.
#' @param X1.variables_BLP A list with three elements controlling the variables that shall be used in the matrix \eqn{X_1} for the BLP regression. The first element of the list, \code{functions_of_Z}, needs to be a subset of \code{c("S", "B", "p")}, where \code{"p"} corresponds to the propensity scores, \code{"B"} to the proxy baseline estimates, and \code{"S"} to the proxy CATE estimates. Default is \code{"B"}. The second element, \code{custom_covariates}, is an optional matrix/data frame of custom covariates that shall be included in \eqn{X_1} (default is \code{NULL}). The third element, \code{fixed_effects}, is a vector of integers that indicates group membership of the observations: For each group, a fixed effect will be added (default is \code{NULL} for no fixed effects). Note that in the final matrix \eqn{X1}, a constant 1 will be silently included so that the regression model has an intercept.
#' @param X1.variables_GATES Same as \code{X1.variables_BLP}, just for the matrix \eqn{X_1} in the GATES regression. Just as in \code{X1.variables_BLP}, a constant 1 will be silently included if no HT transformation is applied so that the GATES regression model has an intercept.
#' @param quantile.cutoffs The cutoff points of quantiles that shall be used for GATES grouping. Default is \code{c(0.25, 0.5, 0.75)}, which corresponds to the four quartiles.
#' @param differences.control_GATES A list with two elements called \code{group.to.subtract.from} and \code{groups.to.be.subtracted}. The first element (\code{group.to.subtract.from}) denotes what shall be the base group to subtract from in the GATES generic targets; either \code{"most"} or \code{"least"}. The second element (\code{groups.to.be.subtracted}) are the groups to be subtracted from \code{group.to.subtract.from}, which is a subset of \eqn{{1,2,...,K}}, where \eqn{K} equals the number of groups. The number of groups should be consistent with the number of groups induced by the argument \code{quantile.cutoffs}.
#' @param differences.control_CLAN Same as \code{differences.control_GATES}, just for the CLAN generic targets.
#' @param vcov.control_BLP A list with two elements called \code{estimator} and \code{arguments}. The element \code{estimator} is a string specifying the covariance matrix estimator to be used in the BLP regression and needs to be a covariance estimator function in the sandwich package (\url{https://cran.r-project.org/web/packages/sandwich/sandwich.pdf}). Recommended estimators are \code{"vcovBS"}, \code{"vcovCL"}, \code{"vcovHAC"}, and \code{"vcovHC"}. Default is \code{"vcovHC"}. The second element \code{arguments} is a list of arguments that shall be passed to the function specified in the first element \code{estimator}. Default leads to the (homoskedastic) ordinary least squares covariance matrix estimator. See the reference manual of the sandwich package for details (\url{https://cran.r-project.org/web/packages/sandwich/vignettes/sandwich.pdf}).
#' @param vcov.control_GATES Same as \code{vcov.control_BLP}, just for GATES regression.
#' @param equal.group.variances_CLAN Logical. If \code{TRUE}, the the two within-group variances of the differences between the CLAN generic targets are assumed to be equal. Default is \code{FALSE}.
#' @param proportion.in.main.set Proportion of samples that shall be in main set. Default is 0.5.
#' @param significance.level Significance level for VEIN. Default is 0.05.
#' @param minimum.variation Minimum variation of the predictions before random noise with distribution \eqn{N(0, var(Y)/20)} is added. Default is \code{1e-05}.
#' @param parallel Logical. If \code{TRUE}, parallel computing will be used. Currently only supported on Unix systems.
#' @param num.cores Number of cores to be used in parallelization (if applicable). Deafult is the number of cores on your machine.
#' @param seed Random seed. Default is \code{NULL} for no random seeding.
#' @param store.learners Logical. If \code{TRUE}, all intermediate results of the learners will be stored. Default is \code{FALSE}. \strong{Warning:} For large data sets and/or many splits in \code{num.splits}, having \code{store.learners = TRUE} might cause memory issues.
#' @param store.splits Logical. If \code{TRUE}, information on the sample splits will be stored. Default is \code{FALSE}.
#'
#' @export
GenericML <- function(Z, D, Y,
                      learner.propensity.score = "constant",
                      learners.genericML,
                      num.splits = 100,
                      Z_CLAN = NULL,
                      HT.transformation          = FALSE,
                      X1.variables_BLP           = list(functions_of_Z = c("B"),
                                                        custom_covariates = NULL,
                                                        fixed_effects = NULL),
                      X1.variables_GATES         = list(functions_of_Z = c("B"),
                                                        custom_covariates = NULL,
                                                        fixed_effects = NULL),
                      quantile.cutoffs           = c(0.25, 0.5, 0.75),
                      differences.control_GATES  = list(group.to.subtract.from = "most",
                                                        groups.to.be.subtracted = 1),
                      differences.control_CLAN   = list(group.to.subtract.from = "most",
                                                        groups.to.be.subtracted = 1),
                      vcov.control_BLP           = list(estimator = "vcovHC",
                                                        arguments = list(type = "const")),
                      vcov.control_GATES         = list(estimator = "vcovHC",
                                                        arguments = list(type = "const")),
                      equal.group.variances_CLAN = FALSE,
                      proportion.in.main.set = 0.5,
                      significance.level = 0.05,
                      minimum.variation = 1e-05,
                      parallel = .Platform$OS.type == "unix",
                      num.cores = parallel::detectCores(),
                      seed = NULL,
                      store.learners = FALSE,
                      store.splits = FALSE){

  ### step 0: input checks ----
  InputChecks_D(D)
  InputChecks_Y(Y)
  InputChecks_Z(Z)
  InputChecks_Z_CLAN(Z_CLAN)
  InputChecks_equal.length3(D, Y, Z)
  InputChecks_X1(X1.variables_BLP)
  InputChecks_X1(X1.variables_GATES)
  InputChecks_vcov.control(vcov.control_BLP)
  InputChecks_vcov.control(vcov.control_GATES)
  InputChecks_differences.control(differences.control_GATES, K = length(quantile.cutoffs) + 1)
  InputChecks_differences.control(differences.control_CLAN, K = length(quantile.cutoffs) + 1)

  if(parallel & .Platform$OS.type != "unix"){
    message("Parallelization is currently only supported on Unix systems (you are using Windows). Therefore, no parallelization will be employed", call. = FALSE)
    parallel <- FALSE

  } # IF

  # render the learners mlr3 environments
  learners <- lapply(1:length(learners.genericML),
                     function(x) get.learner_regr(make.mlr3.environment(learners.genericML[x])))


  ### step 1: compute propensity scores ----
  propensity.scores.obj <- propensity.score_NoChecks(
    Z = Z, D = D, estimator = learner.propensity.score)
  propensity.scores     <- propensity.scores.obj$propensity.scores


  ### step 2: for each ML method, do the generic ML analysis ----

  gen.ml.different.learners <-
    generic.ml.across.learners(Z = Z, D = D, Y = Y,
                               propensity.scores          = propensity.scores,
                               learners                   = learners,
                               learners.names             = learners.genericML,
                               num.splits                 = num.splits,
                               Z_CLAN                     = Z_CLAN,
                               X1.variables_BLP           = X1.variables_BLP,
                               X1.variables_GATES         = X1.variables_GATES,
                               HT.transformation          = HT.transformation,
                               vcov.control_BLP           = vcov.control_BLP,
                               vcov.control_GATES         = vcov.control_GATES,
                               equal.group.variances_CLAN = equal.group.variances_CLAN,
                               proportion.in.main.set     = proportion.in.main.set,
                               quantile.cutoffs           = quantile.cutoffs,
                               differences.control_GATES  = differences.control_GATES,
                               differences.control_CLAN   = differences.control_CLAN,
                               significance.level         = significance.level,
                               minimum.variation          = minimum.variation,
                               parallel                   = parallel,
                               num.cores                  = num.cores,
                               seed                       = seed,
                               store.learners             = store.learners,
                               store.splits               = store.splits)

  # extract the best learners
  best.learners <- get.best.learners(gen.ml.different.learners$generic.targets)


  ### step 3: perform VEIN analysis ----
  vein <- VEIN(gen.ml.different.learners$generic.targets, best.learners)

  # return instance of S3 class 'GenericML'
  return(
    structure(
     list(VEIN = vein,
              best.learners = best.learners,
              propensity.scores = list(estimates = propensity.scores,
                                       mlr3.objects = propensity.scores.obj$mlr3.objects),
              genericML.by.split = gen.ml.different.learners$genericML.by.split,
              splits = gen.ml.different.learners$splits,
              arguments = list(quantile.cutoffs = c(0.25, 0.5, 0.75),
                               proportion.in.main.set = 0.5,
                               significance.level = 0.05,
                               learners.genericML = learners.genericML,
                               learner.propensity.score = learner.propensity.score,
                               num.splits = num.splits,
                               HT.transformation = HT.transformation,
                               X1.variables_BLP = X1.variables_BLP,
                               X1.variables_GATES = X1.variables_GATES)),
     class = "GenericML"))

} # FUN
