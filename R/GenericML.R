#' Generic Machine Learning Inference
#'
#' Performs generic machine learning inference on heterogeneous treatment effects as in \href{https://arxiv.org/abs/1712.04802}{Chernozhukov, Demirer, Duflo and Fernández-Val (2020)} with user-specified machine learning methods. Intended for randomized experiments.
#'
#' @param Z A numeric design matrix that holds the covariates in its columns.
#' @param D A binary vector of treatment assignment. Value one denotes assignment to the treatment group and value zero assignment to the control group.
#' @param Y A numeric vector containing the response variable.
#' @param learners_GenericML A character vector specifying the machine learners to be used for estimating the baseline conditional average (BCA) and conditional average treatment effect (CATE). Either \code{'lasso'}, \code{'random_forest'}, \code{'tree'}, or a custom learner specified with \code{mlr3} syntax. In the latter case, do \emph{not} specify in the \code{mlr3} syntax specification if the learner is a regression learner or classification learner. Example: \code{'mlr3::lrn("ranger", num.trees = 100)'} for a random forest learner with 100 trees. Note that this is a string and the absence of the \code{classif.} or \code{regr.} keywords. See \url{https://mlr3learners.mlr-org.com} for a list of \code{mlr3} learners.
#' @param learner_propensity_score The estimator of the propensity scores. Either a numeric vector (which is then taken as estimates of the propensity scores) or a string specifying the estimator. In the latter case, the string must either be equal to \code{'constant'} (estimates the propensity scores by \code{mean(D)}), \code{'lasso'}, \code{'random_forest'}, \code{'tree'}, or \code{mlr3} syntax. Note that in case of \code{mlr3} syntax, do \emph{not} specify if the learner is a regression learner or classification learner. Example: \code{'mlr3::lrn("ranger", num.trees = 100)'} for a random forest learner with 100 trees. Note that this is a string and the absence of the \code{classif.} or \code{regr.} keywords. See \url{https://mlr3learners.mlr-org.com} for a list of \code{mlr3} learners.
#' @param num_splits Number of sample splits. Default is 100. Must be larger than one. If you want to run \code{GenericML} on a single split, please use \code{\link{GenericML_single}}.
#' @param Z_CLAN A numeric matrix holding variables on which classification analysis (CLAN) shall be performed. CLAN will be performed on each column of the matrix. If \code{NULL} (default), then \code{Z_CLAN = Z}, i.e. CLAN is performed for all variables in \code{Z}.
#' @param HT Logical. If \code{TRUE}, a Horvitz-Thompson (HT) transformation is applied in the BLP and GATES regressions. Default is \code{FALSE}.
#' @param quantile_cutoffs The cutoff points of the quantiles that shall be used for GATES grouping. Default is \code{c(0.25, 0.5, 0.75)}, which corresponds to the four quartiles.
#' @param X1_BLP Specifies the design matrix \eqn{X_1} in the BLP regression. Must be an instance of \code{\link{setup_X1}}. See the documentation of \code{\link{setup_X1}} for details.
#' @param X1_GATES Same as \code{X1_BLP}, just for the GATES regression.
#' @param diff_GATES Specifies the generic targets of GATES. Must be an instance of \code{\link{setup_diff}}. See the documentation of \code{\link{setup_diff}} for details.
#' @param diff_CLAN Same as \code{diff_GATES}, just for the CLAN generic targets.
#' @param vcov_BLP Specifies the covariance matrix estimator in the BLP regression. Must be an instance of \code{\link{setup_vcov}}. See the documentation of \code{\link{setup_vcov}} for details.
#' @param vcov_GATES Same as \code{vcov_BLP}, just for the GATES regression.
#' @param equal_variances_CLAN Logical. If \code{TRUE}, then all within-group variances of the CLAN groups are assumed to be equal. Default is \code{FALSE}. This specification is required for heteroskedasticity-robust variance estimation on the difference of two CLAN generic targets (i.e. variance of the difference of two means). If \code{TRUE} (corresponds to homoskedasticity assumption), the pooled variance is used. If \code{FALSE} (heteroskedasticity), the variance of Welch's t-test is used.
#' @param prop_aux Proportion of samples that shall be in the auxiliary set. Default is 0.5. The number of samples in the auxiliary set will be equal to \code{floor(prop_aux * length(Y))}. If the data set is large, you can save computing time by choosing \code{prop_aux} to be smaller than 0.5.
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
#'   \item{\code{arguments}}{A list of arguments used in the function call.}
#'   }
#'
#' @details
#' The specifications \code{lasso}, \code{random_forest}, and \code{tree} in \code{learners_GenericML} and \code{learner_propensity_score} correspond to the following \code{mlr3} specifications (we omit the keywords \code{classif.} and \code{regr.}). \code{lasso} is a cross-validated Lasso estimator, which corresponds to \code{'mlr3::lrn("cv_glmnet", s = "lambda.min", alpha = 1)'}. \code{random_forest} is a random forest with 500 trees, which corresponds to \code{'mlr3::lrn("ranger", num.trees = 500)'}. \code{tree} is a tree learner, which corresponds to \code{'mlr3::lrn("rpart")'}.
#'
#'
#' @note In an earlier development version, Lucas Kitzmueller alerted us to several minor bugs and proposed fixes. Many thanks to him!
#'
#' @references
#' Chernozhukov V., Demirer M., Duflo E., Fernández-Val I. (2020). \dQuote{Generic Machine Learning Inference on Heterogenous Treatment Effects in Randomized Experiments.} \emph{arXiv preprint arXiv:1712.04802}. URL: \url{https://arxiv.org/abs/1712.04802}.
#'
#' Lang M., Binder M., Richter J., Schratz P., Pfisterer F., Coors S., Au Q., Casalicchio G., Kotthoff L., Bischl B. (2019). \dQuote{mlr3: A Modern Object-Oriented Machine Learning Framework in R.} \emph{Journal of Open Source Software}, \bold{4}(44), 1903. \doi{10.21105/joss.01903}.
#'
#' @seealso
#' \code{\link[=plot.GenericML]{plot}},
#' \code{\link[=print.GenericML]{print}},
#' \code{\link{get_BLP}},
#' \code{\link{get_GATES}},
#' \code{\link{get_CLAN}},
#' \code{\link{setup_X1}},
#' \code{\link{setup_diff}},
#' \code{\link{setup_vcov}},
#' \code{\link{GenericML_single}}
#'
#' @examples
#' if (require("glmnet") && require("ranger")) {
#'
#' ## generate data
#' set.seed(1)
#' n  <- 150                                  # number of observations
#' p  <- 5                                    # number of covariates
#' D  <- rbinom(n, 1, 0.5)                    # random treatment assignment
#' Z  <- matrix(runif(n*p), n, p)             # design matrix
#' Y0 <- as.numeric(Z %*% rexp(p) + rnorm(n)) # potential outcome without treatment
#' Y1 <- 2 + Y0                               # potential outcome under treatment
#' Y  <- ifelse(D == 1, Y1, Y0)               # observed outcome
#'
#' ## column names of Z
#' colnames(Z) <- paste0("V", 1:p)
#'
#' ## specify learners
#' learners <- c("lasso", "mlr3::lrn('ranger', num.trees = 10)")
#'
#' ## glmnet v4.1.3 isn't supported on Solaris, so skip Lasso in this case
#' if(Sys.info()["sysname"] == "SunOS") learners <- learners[-1]
#'
#' ## specify quantile cutoffs (the 4 quartile groups here)
#' quantile_cutoffs <- c(0.25, 0.5, 0.75)
#'
#' ## specify the differenced generic targets of GATES and CLAN
#' # use G4-G1, G4-G2, G4-G3 as differenced generic targets in GATES
#' diff_GATES <- setup_diff(subtract_from = "most",
#'                         subtracted = c(1,2,3))
#' # use G1-G3, G1-G2 as differenced generic targets in CLAN
#' diff_CLAN  <- setup_diff(subtract_from = "least",
#'                          subtracted = c(3,2))
#'
#' ## perform generic ML inference
#' # small number of splits to keep computation time low
#' x <- GenericML(Z, D, Y, learners, num_splits = 2,
#'                quantile_cutoffs = quantile_cutoffs,
#'                diff_GATES = diff_GATES,
#'                diff_CLAN = diff_CLAN,
#'                parallel = FALSE)
#'
#' ## access BLP generic targets for best learner and make plot
#' get_BLP(x, plot = TRUE)
#'
#' ## access GATES generic targets for best learner and make plot
#' get_GATES(x, plot = TRUE)
#'
#' ## access CLAN generic targets for "V1" & best learner and make plot
#' get_CLAN(x, variable = "V1", plot = TRUE)
#'
#' }
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
                      prop_aux                 = 0.5,
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
  InputChecks_num_splits(num_splits)
  stopifnot(is.numeric(quantile_cutoffs))
  stopifnot(0 < min(quantile_cutoffs) & max(quantile_cutoffs) < 1)
  stopifnot(is.logical(equal_variances_CLAN))
  stopifnot(is.logical(HT))
  stopifnot(is.numeric(significance_level) & length(significance_level) == 1)
  stopifnot(0.0 < significance_level & significance_level < 0.5)
  stopifnot(is.numeric(prop_aux) & length(prop_aux) == 1)
  stopifnot(0.0 < prop_aux & prop_aux < 1.0)
  stopifnot(is.numeric(min_variation) & min_variation > 0)
  stopifnot(is.character(learners_GenericML))
  stopifnot(is.character(learner_propensity_score) | is.numeric(learner_propensity_score))

  # if no input provided, set Z_CLAN equal to Z
  if(is.null(Z_CLAN)) Z_CLAN <- Z
  InputChecks_equal.length2(Z, Z_CLAN)

  # set variable names for CLAN
  if(is.null(colnames(Z_CLAN))) colnames(Z_CLAN) <- paste0("V", 1:ncol(Z_CLAN))
  if(any(colnames(Z_CLAN) == "")){

    idx <- which(colnames(Z_CLAN) == "")
    colnames(Z_CLAN)[idx] <- paste0("V", idx)

  } # IF


  if(parallel & !TrueIfUnix()){
    message("Parallelization is currently only supported on Unix systems (you are using Windows). Therefore, no parallelization will be employed", call. = FALSE)
    parallel <- FALSE

  } # IF

  # render the learners mlr3 environments
  learners <- lapply(1:length(learners_GenericML),
                     function(x) get.learner_regr(make.mlr3.environment(learners_GenericML[x])))


  ### step 1: compute propensity scores ----
  if(is.numeric(learner_propensity_score)){

    ## use user-specified propensity scores
    propensity_scores     <- learner_propensity_score
    propensity_scores.obj <- NULL

    # input checks
    InputChecks_equal.length2(propensity_scores, Y)
    if(any(propensity_scores <= 0 | propensity_scores >= 1)){
      stop("User-supplied propensity scores must be contained in the open interval (0,1)",
           call. = FALSE)
    }

  } else{

    ## estimate propensity scores
    stopifnot(is.character(learner_propensity_score))
    stopifnot(length(learner_propensity_score) == 1)

    propensity_scores.obj <-
      propensity_score_NoChecks(
        Z = Z, D = D, estimator = learner_propensity_score)

    propensity_scores <- propensity_scores.obj$estimates

  } # IF

  # check validity of the propensity scores
  InputChecks_propensity_scores(propensity_scores)


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
                               vcov_BLP                   = setup_vcov_align(vcov_BLP),   # align for consistency
                               vcov_GATES                 = setup_vcov_align(vcov_GATES), # align for consistency
                               equal_variances_CLAN       = equal_variances_CLAN,
                               prop_aux                   = prop_aux,
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
          arguments = list(learners_GenericML       = learners_GenericML,
                           learner_propensity_score = learner_propensity_score,
                           num_splits               = num_splits,
                           HT                       = HT,
                           X1_BLP                   = X1_BLP,
                           X1_GATES                 = X1_GATES,
                           quantile_cutoffs         = quantile_cutoffs,
                           diff_GATES               = diff_GATES,
                           diff_CLAN                = diff_CLAN,
                           vcov_BLP                 = vcov_BLP,
                           vcov_GATES               = vcov_GATES,
                           equal_variances_CLAN     = equal_variances_CLAN,
                           prop_aux                 = prop_aux,
                           significance_level       = significance_level,
                           min_variation            = min_variation,
                           parallel                 = parallel,
                           num_cores                = num_cores,
                           seed                     = seed,
                           store_learners           = store_learners,
                           store_splits             = store_splits)),
     class = "GenericML"))

} # FUN
