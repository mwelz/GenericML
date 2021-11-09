#' helper function in case the propensity scores are estimated via mlr3
#'
#' @import mlr3 mlr3learners
#' @noRd
propensity_score_mlr3 <- function(Z, D, learner = "random_forest"){

  # specify the task
  task.propensity_score <- mlr3::TaskClassif$new(id = "propensity_score",
                                                 backend = data.frame(D = as.factor(D), Z),
                                                 target = "D")

  # specify the machine learner
  if(is.environment(learner)){
    learner <- learner
  } else if(learner == "elastic_net"){

    learner <- mlr3::lrn("classif.cv_glmnet", s = "lambda.min")

  } else if(learner == "random_forest"){

    learner <- mlr3::lrn("classif.ranger", num.trees = 500)

  } else if(learner == "tree"){

    learner <- mlr3::lrn("classif.rpart")

  } # END IF


  # specify that the learner predicts Pr(D = 1 | Z)
  learner$predict_type = "prob"

  # fit the learner
  learner$train(task.propensity_score)

  # extract the predictions
  predictions <- learner$predict(task.propensity_score)

  # extract estimations of Pr(D = 1 | Z)
  probs <- predictions$prob

  # return
  return(list(estimates = probs[, colnames(probs) == "1"],
              mlr3_objects = list(task = task.propensity_score,
                                  learner = learner)))

} # END FUN


#' Estimates the propensity scores
#'
#' Estimates the propensity scores \eqn{Pr[D = 1 | Z]} for binary treatment assignment \eqn{D} and covariates \eqn{Z}. Either done by taking the empirical mean of \eqn{D} (which should equal roughly 0.5, since we assume a randomized experiment), or by direct machine learning estimation.
#'
#' @param Z A numeric design matrix that holds the covariates in its columns.
#' @param D A binary vector of treatment assignment. Value one denotes assignment to the treatment group and value zero assignment to the control group.
#' @param estimator The estimator of the propensity scores. Either a numeric vector (which is then taken as estimates of the propensity scores) or a string specifying the estimator. In the latter case, the string must either be equal to \code{'constant'} (estimates the propensity scores by \code{mean(D)}), \code{'elastic_net'}, \code{'random_forest'}, \code{'tree'}, or \code{mlr3} syntax. Note that in case of \code{mlr3} syntax, do \emph{not} specify if the learner is a regression learner or classification learner. Example: \code{'mlr3::lrn("ranger", num.trees = 500)'} for a random forest learner. Note that this is a string and the absence of the \code{classif.} or \code{regr.} keywords. See \url{https://mlr3learners.mlr-org.com} for a list of \code{mlr3} learners.
#'
#' @return
#' An object of class \code{propensity_score}, consisting of the following components:
#' \describe{
#'   \item{\code{estimates}}{A numeric vector of propensity score estimates.}
#'   \item{\code{mlr3_objects}}{\code{mlr3} objects used for estimation. Only non-empty if \code{mlr3} was used.}
#'   }
#'
#' @references
#' Rosenbaum, P.R. and Rubin, D.B. (1983). The Central Role of the Propensity Score in Observational Studies for Causal Effects. \emph{Biometrika}, 70(1):41--55. \url{https://doi.org/10.1093/biomet/70.1.41}.
#'
#' @examples
#' ## generate data
#' library(GenericML)
#' set.seed(1)
#' n  <- 200                        # number of observations
#' p  <- 5                          # number of covariates
#' D  <- rbinom(n, 1, 0.5)          # random treatment assignment
#' Z  <- matrix(runif(n*p), n, p)   # design matrix
#'
#' ## estimate propensity scores
#' propensity_score(Z, D)
#'
#' @export
propensity_score <- function(Z, D, estimator = "constant"){

  # input checks
  InputChecks_D(D)
  InputChecks_Z(Z)
  InputChecks_equal.length2(Z, D)

  # function without input checks
  propensity_score_NoChecks(Z = Z, D = D, estimator = estimator)

} # FUN



#' same as above, but w/o input checks
#'
#' @import mlr3 mlr3learners
#' @noRd
propensity_score_NoChecks <- function(Z, D, estimator = "constant"){

  if(!is.character(estimator)){

    ### case 1: propensity scores are supplied by the user
    if(length(estimator) != length(D)) stop("User-supplied propensity scores in 'estimator' are not of same length as vectors Z and D!")
    if(any(estimator <= 0 | estimator >= 1)) stop("User-supplied propensity scores in 'estimator' must be contained in interval (0,1)!")

    out <- list(estimates = estimator,
                mlr3_objects = NULL)

  } else if(estimator == "constant"){

    ### case 2: propensity scores are estimated by mean of D
    out <- list(estimates = rep(mean(D), length(Z)),
                mlr3_objects = NULL)

  } else{

    ### case 3: propensity scores are estimated by mlr3 (or illegal input)

    # for the following choices of the learner, we require mlr3:
    if(estimator %in% c("elastic_net", "random_forest", "tree") |
       substr(estimator, start = 1, stop = 6) == "mlr3::"){

      out <- propensity_score_mlr3(Z = Z, D = D, learner = make.mlr3.environment(estimator, regr = FALSE))

    } else stop("The argument 'estimator' must be equal to either 'constant', 'elastic_net', random_forest', 'tree', an mlr3 string, or a numeric vector of the same length as Z and D!")

  } # IF

  # return
  return(structure(out, class = "propensity_score"))

} # FUN



#' Performs estimation of the baseline conditional average (BCA), \eqn{E[Y | D=0, Z]}, on the auxiliary sample
#'
#' @param Z A matrix of the covariates.
#' @param D A binary vector of treatment assignment.
#' @param Y The response vector.
#' @param A_set a numerical vector of indices of observations in the auxiliary sample.
#' @param learner A string specifying the machine learner to be used for estimation. Either \code{'elastic_net'}, \code{'random_forest'} (default), or \code{'tree'}. Can alternatively be specified by using \code{mlr3} syntax \emph{without} specification if the learner is a regression learner or classification learner. Example: \code{'mlr3::lrn("ranger", num.trees = 500)'} for a random forest learner. Note that this is a string and the absence of the \code{classif.} or \code{regr.} keywords. See \url{https://mlr3learners.mlr-org.com} for a list of \code{mlr3} learners.
#' @param min_variation Minimum variation of the predictions before random noise with distribution \eqn{N(0, var(Y)/20)} is added. Default is \code{1e-05}.
#'
#' @return An object of the class \code{proxy_BCA}.
#'
#'
#' @export
proxy_BCA <- function(Z, D, Y,
                      A_set,
                      learner = "random_forest",
                      min_variation = 1e-05){

  # input checks
  InputChecks_D(D)
  InputChecks_Y(Y)
  InputChecks_Z(Z)
  InputChecks_equal.length3(D, Y ,Z)

  # call main function
  proxy_BCA_NoChecks(Z = Z, D = D, Y = Y,
                          A_set = A_set,
                          learner = get.learner_regr(learner),
                          min_variation = min_variation)

} # END FUN


#' helper that skips the input checks
#'
#' @import mlr3 mlr3learners
#' @noRd
proxy_BCA_NoChecks <- function(Z, D, Y,
                                    A_set,
                                    learner, # must be mlr3 object
                                    min_variation = 1e-05){

  # specify the task
  task.proxy_BCA.estimator <- mlr3::TaskRegr$new(id = "proxy_BCA",
                                                      backend = data.frame(Y, Z),
                                                      target = "Y")

  # specify that the learner predicts Y
  learner$predict_type = "response"

  # indices of the control units in the auxiliary sample
  A_set.logical <- 1:length(Y) %in% A_set
  idx <- which(A_set.logical & D == 0)

  # fit the learner on the control units in the auxiliary sample
  learner$train(task.proxy_BCA.estimator, row_ids = idx)

  # obtain predictions for Y for all observations
  predictions.obj <- learner$predict(task.proxy_BCA.estimator)
  predictions     <- predictions.obj$response

  # if there is not much variation in the predictions, add Gaussian noise
  if(stats::var(predictions) < min_variation){

    predictions <- predictions +
      stats::rnorm(length(Y), mean = 0, sd = sqrt(stats::var(Y) / 20))

  } # IF

  # return
  return(structure(
    list(estimates  = predictions,
         mlr3_objects = list(task = task.proxy_BCA.estimator,
                             learner = learner)), class = "proxy_BCA"))

} # FUN


#' Performs estimation of the conditional average treatment effect (CATE), \eqn{E[Y | D=1, Z] - E[Y | D=0, Z]}, on the auxiliary sample
#'
#' @param Z A matrix of the covariates.
#' @param D A binary vector of treatment assignment.
#' @param Y The response vector.
#' @param A_set a numerical vector of indices of observations in the auxiliary sample.
#' @param learner A string specifying the machine learner to be used for estimation. Either \code{'elastic_net'}, \code{'random_forest'} (default), or \code{'tree'}. Can alternatively be specified by using \code{mlr3} syntax \emph{without} specification if the learner is a regression learner or classification learner. Example: \code{'mlr3::lrn("ranger", num.trees = 500)'} for a random forest learner. Note that this is a string and the absence of the \code{classif.} or \code{regr.} keywords. See \url{https://mlr3learners.mlr-org.com} for a list of \code{mlr3} learners.
#' @param proxy_BCA A vector of proxy estimates of the baseline estimator BCA, \eqn{E[Y | D=0, Z]}. If \code{NULL}, these will be estimated separately.
#' @param min_variation Minimum variation of the predictions before random noise with distribution \eqn{N(0, var(Y)/20)} is added. Default is \code{1e-05}.
#' @return An object of the class \code{proxy_CATE}.
#'
#' @export
proxy_CATE <- function(Z, D, Y,
                       A_set,
                       learner = "random_forest",
                       proxy_BCA = NULL,
                       min_variation = 1e-05){

  # input checks
  InputChecks_D(D)
  InputChecks_Y(Y)
  InputChecks_Z(Z)
  InputChecks_equal.length3(D, Y ,Z)

  # run main function
  proxy_CATE_NoChecks(Z = Z, D = D, Y = Y,
                      A_set = A_set,
                      learner = get.learner_regr(learner), # must be mlr3 object
                      proxy_BCA = proxy_BCA,
                      min_variation = min_variation)

} # END FUN


#' helper that skips the input checks
#'
#' @import mlr3 mlr3learners
#' @noRd
proxy_CATE_NoChecks <- function(Z, D, Y,
                                A_set,
                                learner = "random_forest",
                                proxy_BCA = NULL,
                                min_variation = 1e-05){

  # indices of the treated units in the auxiliary sample
  A_set.logical <- 1:length(Y) %in% A_set
  idx.auxiliary.treated <- which(A_set.logical & D == 1)

  # specify that the learner predicts Y
  learner$predict_type = "response"

  # specify the task for estimating E[Y | D=1, Z]
  task.proxy_CATE.treated.estimator <- mlr3::TaskRegr$new(id = "cate.treated",
                                                          backend = data.frame(Y, Z),
                                                          target = "Y")

  # prepare the learner
  learner.treated <- learner

  # fit the learner on the treated units in the auxiliary sample
  learner.treated$train(task.proxy_CATE.treated.estimator, row_ids = idx.auxiliary.treated)

  # obtain predictions for E[Y | D=1, Z] for all observations
  predictions.treated.obj <- learner.treated$predict(task.proxy_CATE.treated.estimator)
  predictions.treated     <- predictions.treated.obj$response


  if(!is.null(proxy_BCA)){

    ## if proxy baseline estimates were provided, no need for second estimation
    predictions.controls <- proxy_BCA
    mlr3.controls        <- "Not available as proxy estimates for the baseline were provided"

  } else{

    ## if no proxy baseline estimates were provided, estimate them
    # specify the task for estimating E[Y | D=0, Z]
    task.proxy_CATE.controls.estimator <- mlr3::TaskRegr$new(id = "cate.controls",
                                                             backend = data.frame(Y, Z),
                                                             target = "Y")

    # prepare the learner
    learner.controls <- learner

    # control units in the auxiliary sample
    idx.auxiliary.controls <- which(A_set.logical & D == 0)

    # fit the learner on the control units in the auxiliary sample
    learner.controls$train(task.proxy_CATE.controls.estimator, row_ids = idx.auxiliary.controls)

    # obtain predictions for E[Y | D=0, Z] for all observations
    predictions.controls.obj <- learner.controls$predict(task.proxy_CATE.controls.estimator)
    predictions.controls     <- predictions.controls.obj$response

    # prepare list for output
    mlr3.controls <- list(task = task.proxy_CATE.controls.estimator,
                          learner = learner.controls)

  } # END IF


  # get CATE predictions
  cate.predictions <- predictions.treated - predictions.controls

  # if there is not much variation in the predictions, add Gaussian noise
  if(stats::var(cate.predictions) < min_variation){

    cate.predictions <- cate.predictions +
      stats::rnorm(length(Y), mean = 0, sd = sqrt(stats::var(Y) / 20))

  } # IF


  # prepare output
  estimates <- list(CATE = cate.predictions,
                    Y1 = predictions.treated,
                    Y0 = predictions.controls)

  mlr3_objects <- list(Y1_learner = list(task = task.proxy_CATE.treated.estimator,
                                         learner = learner.treated),
                       Y0_learner = mlr3.controls)

  return(structure(
           list(estimates    = estimates,
                mlr3_objects = mlr3_objects),
         class = "proxy_CATE"))

} # FUN
