#' Generic Machine Learning Inference
#'
#' Performs generic machine learning inference on heterogeneous treatment effects as in \href{https://arxiv.org/abs/1712.04802}{Chernozhukov, Demirer, Duflo and Fernández-Val (2020)} with user-specified machine learning methods. Intended for randomized experiments.
#'
#' @param Z A numeric design matrix that holds the covariates in its columns.
#' @param D A binary vector of treatment assignment. Value one denotes assignment to the treatment group and value zero assignment to the control group.
#' @param Y A numeric vector containing the response variable.
#' @param learners_GenericML A character vector specifying the machine learners to be used for estimating the baseline conditional average (BCA) and conditional average treatment effect (CATE). Either \code{'elastic_net'}, \code{'random_forest'}, \code{'tree'}, or a custom learner specified with \code{mlr3} syntax. In the latter case, do \emph{not} specify in the \code{mlr3} syntax specification if the learner is a regression learner or classification learner. Example: \code{'mlr3::lrn("ranger", num.trees = 500)'} for a random forest learner with 100 trees. Note that this is a string and the absence of the \code{classif.} or \code{regr.} keywords. See \url{https://mlr3learners.mlr-org.com} for a list of \code{mlr3} learners.
#' @param learner_propensity_score The estimator of the propensity scores. Either a numeric vector (which is then taken as estimates of the propensity scores) or a string specifying the estimator. In the latter case, the string must either be equal to \code{'constant'} (estimates the propensity scores by \code{mean(D)}), \code{'elastic_net'}, \code{'random_forest'}, \code{'tree'}, or \code{mlr3} syntax. Note that in case of \code{mlr3} syntax, do \emph{not} specify if the learner is a regression learner or classification learner. Example: \code{'mlr3::lrn("ranger", num.trees = 500)'} for a random forest learner. Note that this is a string and the absence of the \code{classif.} or \code{regr.} keywords. See \url{https://mlr3learners.mlr-org.com} for a list of \code{mlr3} learners.
#' @param num_splits Number of sample splits. Default is 100.
#' @param Z_CLAN A numeric matrix holding variables on which classification analysis (CLAN) shall be performed. CLAN will be performed on each column of the matrix. If \code{NULL} (default), then \code{Z_CLAN = Z}, i.e. CLAN is performed for all variables in \code{Z}.
#' @param HT Logical. If \code{TRUE}, a Horvitz-Thompson (HT) transformation is applied in the BLP and GATES regressions. Default is \code{FALSE}.
#' @param quantile_cutoffs The cutoff points of the quantiles that shall be used for GATES grouping. Default is \code{c(0.25, 0.5, 0.75)}, which corresponds to the four quartiles.
#' @param X1_BLP Specifies the design matrix \eqn{X_1} in the BLP regression. See the documentation of \code{\link{setup_X1}} for details.
#' @param X1_GATES Same as \code{X1_BLP}, just for the GATES regression.
#' @param diff_GATES Specifies the generic targets of GATES. See the documentation of \code{\link{setup_diff}} for details.
#' @param diff_CLAN Same as \code{diff_GATES}, just for the CLAN generic targets.
#' @param vcov_BLP Specifies the covariance matrix estimator in the BLP regression. See the documentation of \code{\link{setup_vcov}} for details.
#' @param vcov_GATES Same as \code{vcov_BLP}, just for the GATES regression.
#' @param equal_variances_CLAN Logical. If \code{TRUE}, then all within-group variances of the CLAN groups are assumed to be equal. Default is \code{FALSE}. This specification is required for heteroskedasticity-robust variance estimation on the difference of two CLAN generic targets (i.e. variance of the difference of two means). If \code{TRUE} (corresponds to homoskedasticity assumption), the pooled variance is used. If \code{FALSE} (heteroskedasticity), the variance of Welch's t-test is used.
#' @param prop_main Proportion of samples that shall be in main set. Default is 0.5.
#' @param significance_level Significance level for VEIN. Default is 0.05.
#' @param min_variation Specifies a threshold for the minimum variation of the BCA/CATE predictions. If the variation of a BCA/CATE prediction falls below this threshold, random noise with distribution \eqn{N(0, var(Y)/20)} is added to it. Default is \code{1e-05}.
#' @param parallel Logical. If \code{TRUE}, parallel computing will be used. Currently only supported for Unix systems.
#' @param num_cores Number of cores to be used in parallelization (if applicable). Default is the number of cores of the user's machine.
#' @param seed Random seed. Default is \code{NULL} for no random seeding.
#' @param store_learners Logical. If \code{TRUE}, all intermediate results of the learners will be stored. That is, for each learner and each split, all BCA and CATE predictions as well as all BLP, GATES, CLAN, and \eqn{\Lambda} estimates will be stored. Default is \code{FALSE}.
#' @param store_splits Logical. If \code{TRUE} (default), the sample splits will be stored.
#'
#' @return
#' An object of class \code{GenericML}. On this object, we recommend to use the accessor functions \code{\link{get_BLP}}, \code{\link{get_GATES}}, and \code{\link{get_CLAN}} to extract the results of the analyses of BLP, GATES, and CLAN, respectively. An object of class \code{GenericML} contains the following components:
#' \describe{
#'   \item{\code{VEIN}}{A list containing two sub-lists called \code{best_learners} and \code{all_learners}, respectively. Each of these two sub-lists contains the inferential VEIN results on the generic targets of the BLP, GATES, and CLAN analyses. \code{all_learners} does this for all learners specified in the argument \code{learners_GenericML}, \code{best_learners} only for the corresponding best learners. Which learner is best for which analysis is assessed by the \eqn{\Lambda} criteria discussed in Sections 5.2 and 5.3 of the paper.}
#'   \item{\code{best}}{A list containing information on the evaluation of which learner is the best for which analysis. Contains four components. The first three contain the name of the best learner for BLP, GATES, and CLAN, respectively. The fourth component, \code{overview}, contains the two \eqn{\Lambda} criteria used to determine the best learners (discussed in Sections 5.2 and 5.3 of the paper).}
#'   \item{\code{propensity_scores}}{The propensity score estimates as well as the \code{mlr3} objects used to estimate them (if \code{mlr3} was used for estimation).}
#'   \item{\code{GenericML_single}}{Only nonempty if \code{store_learners = TRUE}. Contains all intermediate results of each learners for each split. That is, for a given learner (first level of the list) and split (second level),  objects of classes \code{\link{BLP}}, \code{\link{GATES}}, \code{\link{CLAN}}, \code{\link{proxy_BCA}}, \code{\link{proxy_CATE}} as well as the \eqn{\Lambda} criteria (\code{"best"})) are listed, which were computed with the given learner and split.}
#'   \item{\code{splits}}{Only nonempty if \code{store_splits = TRUE}. Contains a character matrix of dimension \code{length(Y)} by \code{num_splits}. Contains the group membership (main or auxiliary) of each observation (rows) in each split (columns). \code{"M"} denotes the main set, \code{"A"} the auxiliary set.}
#'   \item{\code{arguments}}{A list of all arguments used in the function call.}
#'   }
#'
#' @note In an earlier development version, Lucas Kitzmueller alerted us to several minor bugs and proposed fixes. Many thanks to him!
#'
#' @references
#' Chernozhukov, V., Demirer, M., Duflo, E., and Fernández-Val, I. (2021). Generic Machine Learning Inference on Heterogenous Treatment Effects in Randomized Experiments. \href{https://arxiv.org/abs/1712.04802}{\emph{arXiv preprint arXiv:1712.04802}}.
#'
#' @seealso \code{\link[=plot.GenericML]{plot}},
#' \code{\link[=print.GenericML]{print}},
#' \code{\link{get_BLP}},
#' \code{\link{get_GATES}},
#' \code{\link{get_CLAN}},
#' \code{\link{setup_X1}},
#' \code{\link{setup_diff}},
#' \code{\link{setup_vcov}},
#' \code{\link{GenericML_single}}
#'
#'
#' @export
GenericML <- function(Z, D, Y,
                      learners_GenericML,
                      learner_propensity_score = "constant",
                      num_splits               = 100,
                      Z_CLAN                   = NULL,
                      HT                       = FALSE,
                      quantile_cutoffs         = c(0.25, 0.5, 0.75),
                      X1_BLP                   = setup_X1(),
                      X1_GATES                 = setup_X1(),
                      diff_GATES               = setup_diff(),
                      diff_CLAN                = setup_diff(),
                      vcov_BLP                 = setup_vcov(),
                      vcov_GATES               = setup_vcov(),
                      equal_variances_CLAN     = FALSE,
                      prop_main                = 0.5,
                      significance_level       = 0.05,
                      min_variation            = 1e-05,
                      parallel                 = TrueIfUnix(),
                      num_cores                = parallel::detectCores(),
                      seed                     = NULL,
                      store_learners           = FALSE,
                      store_splits             = TRUE){

  ### step 0: input checks ----
  InputChecks_D(D)
  InputChecks_Y(Y)
  InputChecks_Z(Z)
  InputChecks_Z_CLAN(Z_CLAN)
  InputChecks_equal.length3(D, Y, Z)
  InputChecks_X1(X1_BLP, length(Y))
  InputChecks_X1(X1_GATES, length(Y))
  InputChecks_vcov.control(vcov_BLP)
  InputChecks_vcov.control(vcov_GATES)
  InputChecks_diff(diff_GATES, K = length(quantile_cutoffs) + 1)
  InputChecks_diff(diff_CLAN, K = length(quantile_cutoffs) + 1)
  stopifnot(is.numeric(quantile_cutoffs))
  stopifnot(is.logical(equal_variances_CLAN))
  stopifnot(is.logical(HT))
  stopifnot(is.numeric(significance_level))
  stopifnot(is.numeric(prop_main))

  if(parallel & !TrueIfUnix()){
    message("Parallelization is currently only supported on Unix systems (you are using Windows). Therefore, no parallelization will be employed", call. = FALSE)
    parallel <- FALSE

  } # IF

  # render the learners mlr3 environments
  learners <- lapply(1:length(learners_GenericML),
                     function(x) get.learner_regr(make.mlr3.environment(learners_GenericML[x])))


  ### step 1: compute propensity scores ----
  propensity_scores.obj <- propensity_score_NoChecks(
    Z = Z, D = D, estimator = learner_propensity_score)
  propensity_scores     <- propensity_scores.obj$estimates


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
  best.learners <- get.best.learners(gen.ml.different.learners$generic_targets)


  ### step 3: perform VEIN analysis ----
  vein <- VEIN(gen.ml.different.learners$generic_targets, best.learners)

  # return instance of S3 class 'GenericML'
  return(
    structure(
     list(VEIN = vein,
          best = best.learners,
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
