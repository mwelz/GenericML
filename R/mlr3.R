#' Propensity score estimation
#'
#' Estimates the propensity scores \eqn{Pr[D = 1 | Z]} for binary treatment assignment \eqn{D} and covariates \eqn{Z}. Either done by taking the empirical mean of \eqn{D} (which should equal roughly 0.5, since we assume a randomized experiment), or by direct machine learning estimation.
#'
#' @param Z A numeric design matrix that holds the covariates in its columns.
#' @param D A binary vector of treatment assignment. Value one denotes assignment to the treatment group and value zero assignment to the control group.
#' @param estimator Character specifying the estimator. Must either be equal to \code{'constant'} (estimates the propensity scores by \code{mean(D)}), \code{'lasso'}, \code{'random_forest'}, \code{'tree'}, or \code{mlr3} syntax. Note that in case of \code{mlr3} syntax, do \emph{not} specify if the learner is a regression learner or classification learner. Example: \code{'mlr3::lrn("ranger", num.trees = 500)'} for a random forest learner. Note that this is a string and the absence of the \code{classif.} or \code{regr.} keywords. See \url{https://mlr3learners.mlr-org.com} for a list of \code{mlr3} learners.
#'
#' @return
#' An object of class \code{propensity_score}, consisting of the following components:
#' \describe{
#'   \item{\code{estimates}}{A numeric vector of propensity score estimates.}
#'   \item{\code{mlr3_objects}}{\code{mlr3} objects used for estimation. Only non-empty if \code{mlr3} was used.}
#'   }
#'
#' @details
#' The specifications \code{lasso}, \code{random_forest}, and \code{tree} in \code{estimator} correspond to the following \code{mlr3} specifications (we omit the keywords \code{classif.} and \code{regr.}). \code{lasso} is a cross-validated Lasso estimator, which corresponds to \code{'mlr3::lrn("cv_glmnet", s = "lambda.min", alpha = 1)'}. \code{random_forest} is a random forest with 500 trees, which corresponds to \code{'mlr3::lrn("ranger", num.trees = 500)'}. \code{tree} is a tree learner, which corresponds to \code{'mlr3::lrn("rpart")'}.
#'
#' @references
#' Rosenbaum P.R., Rubin D.B. (1983). \dQuote{The Central Role of the Propensity Score in Observational Studies for Causal Effects.} \emph{Biometrika}, \bold{70}(1), 41--55. \doi{10.1093/biomet/70.1.41}.
#'
#' Lang M., Binder M., Richter J., Schratz P., Pfisterer F., Coors S., Au Q., Casalicchio G., Kotthoff L., Bischl B. (2019). \dQuote{mlr3: A Modern Object-Oriented Machine Learning Framework in R.} \emph{Journal of Open Source Software}, \bold{4}(44), 1903. \doi{10.21105/joss.01903}.
#'
#' @examples
#' ## generate data
#' set.seed(1)
#' n  <- 100                        # number of observations
#' p  <- 5                          # number of covariates
#' D  <- rbinom(n, 1, 0.5)          # random treatment assignment
#' Z  <- matrix(runif(n*p), n, p)   # design matrix
#'
#' ## estimate propensity scores via mean(D)...
#' propensity_score(Z, D, estimator = "constant")
#'
#' ## ... and via SVM with cache size 40
#' if(require("e1071")){
#'   propensity_score(Z, D,
#'    estimator = 'mlr3::lrn("svm", cachesize = 40)')
#' }
#'
#' @export
propensity_score <- function(Z, D, estimator = "constant"){

  # input checks
  InputChecks_D(D)
  InputChecks_Z(Z)
  InputChecks_equal.length2(Z, D)
  stopifnot(is.character(estimator))
  stopifnot(length(estimator) == 1)

  # function without input checks
  propensity_score_NoChecks(Z = Z, D = D, estimator = estimator)

} # FUN



#' propensity_score(), but w/o input checks
#'
#' @import mlr3 mlr3learners
#' @noRd
propensity_score_NoChecks <- function(Z, D, estimator = "constant"){

  if(estimator == "constant"){

    out <- list(estimates = rep(mean(D), length(D)),
                mlr3_objects = NULL)

  } else if(estimator %in% c("lasso", "random_forest", "tree") |
            substr(estimator, start = 1, stop = 6) == "mlr3::"){

    out <- propensity_score_mlr3(Z = Z, D = D,
                                 learner = make.mlr3.environment(estimator, regr = FALSE))

  } else stop(paste0("The argument 'estimator' must be equal to either",
              "'constant', 'lasso', random_forest', 'tree', or an mlr3 string"))

  # return
  return(structure(out, class = "propensity_score"))

} # FUN


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
  } else if(learner == "lasso"){

    learner <- mlr3::lrn("classif.cv_glmnet", s = "lambda.min", alpha = 1)

  } else if(learner == "random_forest"){

    learner <- mlr3::lrn("classif.ranger", num.trees = 500)

  } else if(learner == "tree"){

    learner <- mlr3::lrn("classif.rpart")

  } # END IF


  # specify that the learner predicts Pr(D = 1 | Z)
  learner$predict_type = "prob"

  # ensure that mlr3 only uses one thread (default in mlr3)
  invisible(mlr3::set_threads(x = learner, n = 1L))

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



#' Baseline Conditional Average
#'
#' Proxy estimation of the Baseline Conditional Average (BCA), defined by \eqn{E[Y | D=0, Z]}. Estimation is done on the auxiliary sample, but BCA predictions are made for all observations.
#'
#' @param Z A numeric design matrix that holds the covariates in its columns.
#' @param D A binary vector of treatment assignment. Value one denotes assignment to the treatment group and value zero assignment to the control group.
#' @param Y A numeric vector containing the response variable.
#' @param A_set A numerical vector of the indices of the observations in the auxiliary sample.
#' @param learner A string specifying the machine learner for the estimation. Either \code{'lasso'}, \code{'random_forest'}, \code{'tree'}, or a custom learner specified with \code{mlr3} syntax. In the latter case, do \emph{not} specify in the \code{mlr3} syntax specification if the learner is a regression learner or classification learner. Example: \code{'mlr3::lrn("ranger", num.trees = 100)'} for a random forest learner with 100 trees. Note that this is a string and the absence of the \code{classif.} or \code{regr.} keywords. See \url{https://mlr3learners.mlr-org.com} for a list of \code{mlr3} learners.
#' @param min_variation Specifies a threshold for the minimum variation of the predictions. If the variation of a BCA prediction falls below this threshold, random noise with distribution \eqn{N(0, var(Y)/20)} is added to it. Default is \code{1e-05}.
#'
#' @return
#' An object of class \code{proxy_BCA}, consisting of the following components:
#' \describe{
#'   \item{\code{estimates}}{A numeric vector of BCA estimates of each observation.}
#'   \item{\code{mlr3_objects}}{\code{mlr3} objects used for estimation.}
#'   }
#'
#' @details
#' The specifications \code{lasso}, \code{random_forest}, and \code{tree} in \code{learner} correspond to the following \code{mlr3} specifications (we omit the keywords \code{classif.} and \code{regr.}). \code{lasso} is a cross-validated Lasso estimator, which corresponds to \code{'mlr3::lrn("cv_glmnet", s = "lambda.min", alpha = 1)'}. \code{random_forest} is a random forest with 500 trees, which corresponds to \code{'mlr3::lrn("ranger", num.trees = 500)'}. \code{tree} is a tree learner, which corresponds to \code{'mlr3::lrn("rpart")'}.
#'
#' @references
#' Chernozhukov V., Demirer M., Duflo E., Fernández-Val I. (2020). \dQuote{Generic Machine Learning Inference on Heterogenous Treatment Effects in Randomized Experiments.} \emph{arXiv preprint arXiv:1712.04802}. URL: \url{https://arxiv.org/abs/1712.04802}.
#'
#' Lang M., Binder M., Richter J., Schratz P., Pfisterer F., Coors S., Au Q., Casalicchio G., Kotthoff L., Bischl B. (2019). \dQuote{mlr3: A Modern Object-Oriented Machine Learning Framework in R.} \emph{Journal of Open Source Software}, \bold{4}(44), 1903. \doi{10.21105/joss.01903}.
#'
#' @seealso \code{\link{proxy_CATE}}
#'
#' @examples
#' if(require("ranger")){
#' ## generate data
#' set.seed(1)
#' n  <- 150                                  # number of observations
#' p  <- 5                                    # number of covariates
#' D  <- rbinom(n, 1, 0.5)                    # random treatment assignment
#' Z  <- matrix(runif(n*p), n, p)             # design matrix
#' Y0 <- as.numeric(Z %*% rexp(p) + rnorm(n)) # potential outcome without treatment
#' Y1 <- 2 + Y0                               # potential outcome under treatment
#' Y  <- ifelse(D == 1, Y1, Y0)               # observed outcome
#' A_set <- sample(1:n, size = n/2)           # auxiliary set
#'
#' ## BCA predictions via random forest
#' proxy_BCA(Z, D, Y, A_set, learner = "mlr3::lrn('ranger', num.trees = 10)")
#' }
#'
#' @export
proxy_BCA <- function(Z, D, Y,
                      A_set,
                      learner,
                      min_variation = 1e-05){

  # input checks
  InputChecks_D(D)
  InputChecks_Y(Y)
  InputChecks_Z(Z)
  InputChecks_equal.length3(D, Y ,Z)
  stopifnot(is.character(learner))
  stopifnot(length(learner) == 1)
  InputChecks_index_set(A_set, nrow(Z))
  stopifnot(is.numeric(min_variation) & min_variation > 0)


  # call main function
  proxy_BCA_NoChecks(Z = Z, D = D, Y = Y,
                          A_set = A_set,
                          learner = get.learner_regr(make.mlr3.environment(learner)), # must be mlr3 object
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

  # ensure that mlr3 only uses one thread (default in mlr3)
  invisible(mlr3::set_threads(x = learner, n = 1L))

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


#' Conditional Average Treatment Effect
#'
#' Proxy estimation of the Conditional Average Treatment Effect (CATE), defined by \eqn{E[Y | D=1, Z] - E[Y | D=0, Z]}. Estimation is done on the auxiliary sample, but CATE predictions are made for all observations.
#'
#' @param Z A numeric design matrix that holds the covariates in its columns.
#' @param D A binary vector of treatment assignment. Value one denotes assignment to the treatment group and value zero assignment to the control group.
#' @param Y A numeric vector containing the response variable.
#' @param A_set A numerical vector of the indices of the observations in the auxiliary sample.
#' @param learner A string specifying the machine learner for the estimation. Either \code{'lasso'}, \code{'random_forest'}, \code{'tree'}, or a custom learner specified with \code{mlr3} syntax. In the latter case, do \emph{not} specify in the \code{mlr3} syntax specification if the learner is a regression learner or classification learner. Example: \code{'mlr3::lrn("ranger", num.trees = 100)'} for a random forest learner with 100 trees. Note that this is a string and the absence of the \code{classif.} or \code{regr.} keywords. See \url{https://mlr3learners.mlr-org.com} for a list of \code{mlr3} learners.
#' @param proxy_BCA A vector of proxy estimates of the baseline conditional average, BCA, \eqn{E[Y | D=0, Z]}. If \code{NULL}, these will be estimated separately.
#' @param min_variation Minimum variation of the predictions before random noise with distribution \eqn{N(0, var(Y)/20)} is added. Default is \code{1e-05}.
#'
#' @return
#' An object of class \code{proxy_CATE}, consisting of the following components:
#' \describe{
#'   \item{\code{estimates}}{A numeric vector of CATE estimates of each observation.}
#'   \item{\code{mlr3_objects}}{\code{mlr3} objects used for estimation of \eqn{E[Y | D=1, Z]} (\code{Y1_learner}) and \eqn{E[Y | D=0, Z]} (\code{Y0_learner}). The latter is not available if \code{proxy_BCA = NULL}.}
#'   }
#'
#' @details
#' The specifications \code{lasso}, \code{random_forest}, and \code{tree} in \code{learner} correspond to the following \code{mlr3} specifications (we omit the keywords \code{classif.} and \code{regr.}). \code{lasso} is a cross-validated Lasso estimator, which corresponds to \code{'mlr3::lrn("cv_glmnet", s = "lambda.min", alpha = 1)'}. \code{random_forest} is a random forest with 500 trees, which corresponds to \code{'mlr3::lrn("ranger", num.trees = 500)'}. \code{tree} is a tree learner, which corresponds to \code{'mlr3::lrn("rpart")'}.
#'
#' @references
#' Chernozhukov V., Demirer M., Duflo E., Fernández-Val I. (2020). \dQuote{Generic Machine Learning Inference on Heterogenous Treatment Effects in Randomized Experiments.} \emph{arXiv preprint arXiv:1712.04802}. URL: \url{https://arxiv.org/abs/1712.04802}.
#'
#' Lang M., Binder M., Richter J., Schratz P., Pfisterer F., Coors S., Au Q., Casalicchio G., Kotthoff L., Bischl B. (2019). \dQuote{mlr3: A Modern Object-Oriented Machine Learning Framework in R.} \emph{Journal of Open Source Software}, \bold{4}(44), 1903. \doi{10.21105/joss.01903}.
#'
#' @seealso \code{\link{proxy_BCA}}
#'
#' @examples
#' if(require("ranger")){
#' ## generate data
#' set.seed(1)
#' n  <- 150                                  # number of observations
#' p  <- 5                                    # number of covariates
#' D  <- rbinom(n, 1, 0.5)                    # random treatment assignment
#' Z  <- matrix(runif(n*p), n, p)             # design matrix
#' Y0 <- as.numeric(Z %*% rexp(p) + rnorm(n)) # potential outcome without treatment
#' Y1 <- 2 + Y0                               # potential outcome under treatment
#' Y  <- ifelse(D == 1, Y1, Y0)               # observed outcome
#' A_set <- sample(1:n, size = n/2)           # auxiliary set
#'
#' ## CATE predictions via random forest
#' proxy_CATE(Z, D, Y, A_set, learner = "mlr3::lrn('ranger', num.trees = 10)")
#' }
#'
#' @export
proxy_CATE <- function(Z, D, Y,
                       A_set,
                       learner,
                       proxy_BCA = NULL,
                       min_variation = 1e-05){


  # input checks
  InputChecks_D(D)
  InputChecks_Y(Y)
  InputChecks_Z(Z)
  InputChecks_equal.length3(D, Y ,Z)
  stopifnot(is.character(learner))
  stopifnot(length(learner) == 1)
  InputChecks_index_set(A_set, nrow(Z))
  stopifnot(is.numeric(min_variation) & min_variation > 0)

  if(!is.null(proxy_BCA)){
    stopifnot(is.numeric(proxy_BCA) & is.vector(proxy_BCA))
    InputChecks_equal.length2(proxy_BCA, Y)
  } # IF

  # run main function
  proxy_CATE_NoChecks(Z = Z, D = D, Y = Y,
                      A_set = A_set,
                      learner = get.learner_regr(make.mlr3.environment(learner)), # must be mlr3 object
                      proxy_BCA = proxy_BCA,
                      min_variation = min_variation)

} # END FUN


#' helper that skips the input checks
#'
#' @import mlr3 mlr3learners
#' @noRd
proxy_CATE_NoChecks <- function(Z, D, Y,
                                A_set,
                                learner, # must be mlr3 object
                                proxy_BCA = NULL,
                                min_variation = 1e-05){

  # indices of the treated units in the auxiliary sample
  A_set.logical <- 1:length(Y) %in% A_set
  idx.auxiliary.treated <- which(A_set.logical & D == 1)

  # specify that the learner predicts Y
  learner$predict_type = "response"

  # ensure that mlr3 only uses one thread (default in mlr3)
  invisible(mlr3::set_threads(x = learner, n = 1L))

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
