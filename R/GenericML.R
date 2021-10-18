#' Generic Machine Learning Inference
#'
#' Performs generic machine learning inference as in Chernozhukov, Demirer, Duflo and Fernández-Val (2020). Link to working paper: \url{https://arxiv.org/abs/1712.04802}.
#'
#' @param Z A matrix or data frame of the covariates.
#' @param D A binary vector of treatment assignment.
#' @param Y The response vector.
#' @param learners_GenericML A vector of strings specifying the machine learners to be used for estimating the BCA and CATE. Either \code{'elastic.net'}, \code{'random.forest'}, or \code{'tree'}. Can alternatively be specified by using \code{mlr3} syntax \emph{without} specification if the learner is a regression learner or classification learner. Example: \code{'mlr3::lrn("ranger", num.trees = 500)'} for a random forest learner. Note that this is a string and the absence of the \code{classif.} or \code{regr.} keyords. See \url{https://mlr3learners.mlr-org.com} for a list of \code{mlr3} learners.
#' @param learner_propensity_score The estimator for the propensity scores. Either a numeric vector (which is then taken as estimates of the propensity scores) or a string specifying the estimator. The string must either be equal to \code{'constant'} (estimates the propensity scores by \code{mean(D)}), \code{'elastic.net'}, \code{'random.forest'}, \code{'tree'}, or \code{mlr3} syntax. Note that in case of \code{mlr3} syntax, do \emph{not} specify if the learner is a regression learner or classification learner; Example: \code{'mlr3::lrn("ranger", num.trees = 500)'} for a random forest learner. Note that this is a string and the absence of the \code{classif.} or \code{regr.} keyords.
#' @param num_splits Number of sample splits. Default is 100.
#' @param Z_CLAN A matrix of variables that shall be considered for the CLAN. If \code{NULL} (default), then \code{Z_CLAN = Z}, i.e. CLAN is performed for all variables in \code{Z}.
#' @param HT Logical. If TRUE, a HT transformation is applied in BLP and GATES. Default is FALSE.
#' @param X1_BLP Specifies the design matrix \eqn{X_1} in the BLP regression. See the documentation of \code{\link{setup_X1}} for details.
#' @param X1_GATES Same as \code{X1_BLP}, just for the the GATES regression.
#' @param quantile_cutoffs The cutoff points of quantiles that shall be used for GATES grouping. Default is \code{c(0.25, 0.5, 0.75)}, which corresponds to the four quartiles.
#' @param diff_GATES Specifies the generic targets of GATES. See the documentation of \code{\link{setup_diff}} for details.
#' @param diff_CLAN Same as \code{diff_GATES}, just for the CLAN generic targets.
#' @param vcov_BLP Specifies the covariance matrix estimator in the BLP regression. See the documentation of \code{\link{setup_vcov}} for details.
#' @param vcov_GATES Same as \code{vcov_BLP}, just for the GATES regression.
#' @param equal_variances_CLAN Logical. If \code{TRUE}, the the two within-group variances of the differences between the CLAN generic targets are assumed to be equal. Default is \code{FALSE}.
#' @param prop_main Proportion of samples that shall be in main set. Default is 0.5.
#' @param significance_level Significance level for VEIN. Default is 0.05.
#' @param min_variation Minimum variation of the predictions before random noise with distribution \eqn{N(0, var(Y)/20)} is added. Default is \code{1e-05}.
#' @param parallel Logical. If \code{TRUE}, parallel computing will be used. Currently only supported on Unix systems.
#' @param num_cores Number of cores to be used in parallelization (if applicable). Deafult is the number of cores on your machine.
#' @param seed Random seed. Default is \code{NULL} for no random seeding.
#' @param store_learners Logical. If \code{TRUE}, all intermediate results of the learners will be stored. Default is \code{FALSE}. \strong{Warning:} For large data sets and/or many splits in \code{num_splits}, having \code{store_learners = TRUE} might cause memory issues.
#' @param store_splits Logical. If \code{TRUE}, information on the sample splits will be stored. Default is \code{FALSE}.
#'
#' @export
GenericML <- function(Z, D, Y,
                      learners_GenericML,
                      learner_propensity_score = "constant",
                      num_splits               = 100,
                      Z_CLAN                   = NULL,
                      HT                       = FALSE,
                      X1_BLP                   = setup_X1(),
                      X1_GATES                 = setup_X1(),
                      quantile_cutoffs         = c(0.25, 0.5, 0.75),
                      diff_GATES               = setup_diff(),
                      diff_CLAN                = setup_diff(),
                      vcov_BLP                 = setup_vcov(),
                      vcov_GATES               = setup_vcov(),
                      equal_variances_CLAN     = FALSE,
                      prop_main                = 0.5,
                      significance_level       = 0.05,
                      min_variation            = 1e-05,
                      parallel                 = .Platform$OS.type == "unix",
                      num_cores                = parallel::detectCores(),
                      seed                     = NULL,
                      store_learners           = FALSE,
                      store_splits             = FALSE){

  ### step 0: input checks ----
  InputChecks_D(D)
  InputChecks_Y(Y)
  InputChecks_Z(Z)
  InputChecks_Z_CLAN(Z_CLAN)
  InputChecks_equal.length3(D, Y, Z)
  InputChecks_X1(X1_BLP)
  InputChecks_X1(X1_GATES)
  InputChecks_vcov.control(vcov_BLP)
  InputChecks_vcov.control(vcov_GATES)
  InputChecks_diff(diff_GATES, K = length(quantile_cutoffs) + 1)
  InputChecks_diff(diff_CLAN, K = length(quantile_cutoffs) + 1)

  if(parallel & .Platform$OS.type != "unix"){
    message("Parallelization is currently only supported on Unix systems (you are using Windows). Therefore, no parallelization will be employed", call. = FALSE)
    parallel <- FALSE

  } # IF

  # render the learners mlr3 environments
  learners <- lapply(1:length(learners_GenericML),
                     function(x) get.learner_regr(make.mlr3.environment(learners_GenericML[x])))


  ### step 1: compute propensity scores ----
  propensity_scores.obj <- propensity_score_NoChecks(
    Z = Z, D = D, estimator = learner_propensity_score)
  propensity_scores     <- propensity_scores.obj$propensity_scores


  ### step 2: for each ML method, do the generic ML analysis ----

  gen.ml.different.learners <-
    generic.ml.across.learners(Z = Z, D = D, Y = Y,
                               propensity_scores          = propensity_scores,
                               learners                   = learners,
                               learners.names             = learners_GenericML,
                               num_splits                 = num_splits,
                               Z_CLAN                     = Z_CLAN,
                               X1_BLP                     = X1_BLP,
                               X1_GATES                   = X1_GATES,
                               HT                         = HT,
                               vcov_BLP                   = vcov_BLP,
                               vcov_GATES                 = vcov_GATES,
                               equal_variances_CLAN       = equal_variances_CLAN,
                               prop_main                  = prop_main,
                               quantile_cutoffs           = quantile_cutoffs,
                               diff_GATES                 = diff_GATES,
                               diff_CLAN                  = diff_CLAN,
                               significance_level         = significance_level,
                               min_variation              = min_variation,
                               parallel                   = parallel,
                               num_cores                  = num_cores,
                               seed                       = seed,
                               store_learners             = store_learners,
                               store_splits               = store_splits)

  # extract the best learners
  best.learners <- get.best.learners(gen.ml.different.learners$generic.targets)


  ### step 3: perform VEIN analysis ----
  vein <- VEIN(gen.ml.different.learners$generic.targets, best.learners)

  # return instance of S3 class 'GenericML'
  return(
    structure(
     list(VEIN = vein,
          best_learners = best.learners,
          propensity_scores = list(estimates = propensity_scores,
                                   mlr3_objects = propensity_scores.obj$mlr3_objects),
          GenericML_single = gen.ml.different.learners$genericML.by.split,
          splits = gen.ml.different.learners$splits,
          arguments = list(Z = Z, D = D, Y = Y,
                           learners_GenericML = learners_GenericML,
                           learner_propensity_score = learner_propensity_score,
                           num_splits               = num_splits,
                           Z_CLAN                   = Z_CLAN,
                           HT                       = HT,
                           X1_BLP                   = X1_BLP,
                           X1_GATES                 = X1_GATES,
                           quantile_cutoffs         = quantile_cutoffs,
                           diff_GATES               = diff_GATES,
                           diff_CLAN                = diff_CLAN,
                           vcov_BLP                 = vcov_BLP,
                           vcov_GATES               = vcov_GATES,
                           equal_variances_CLAN     = equal_variances_CLAN,
                           prop_main                = prop_main,
                           significance_level       = significance_level,
                           min_variation            = min_variation,
                           parallel                 = parallel,
                           num_cores                = num_cores,
                           seed                     = seed,
                           store_learners           = store_learners,
                           store_splits             = store_splits)),
     class = "GenericML"))

} # FUN
