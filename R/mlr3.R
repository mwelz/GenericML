#' helper function in case the propensity scores are estimated via mlr3
#'
#' @import mlr3 mlr3learners
#' @noRd
propensity_score_mlr3 <- function(Z, D, learner = "random.forest"){

  # specify the task
  task.propensity_score <- mlr3::TaskClassif$new(id = "propensity_score",
                                                 backend = data.frame(D = as.factor(D), Z),
                                                 target = "D")

  # specify the machine learner
  if(is.environment(learner)){
    learner <- learner
  } else if(learner == "elastic.net"){

    learner <- mlr3::lrn("classif.cv_glmnet", s = "lambda.min")

  } else if(learner == "random.forest"){

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
  return(list(propensity_scores = probs[, colnames(probs) == "1"],
              mlr3_objects = list(task = task.propensity_score,
                                  learner = learner)))

} # END FUN


#' Estimates the propensity score
#'
#' @param Z a matrix or data frame of covariates
#' @param D a binary vector of treatment status
#' @param estimator the estimator to be used. Either a numeric vector (which is then taken as estimates of the propensity scores) or a string specifying the estimator. The string must either be equal to 'constant' (estimates the propensity scores by mean(D)), 'elastic.net', 'random.forest', 'tree', or mlr3 syntax. Example for the latter: mlr3::lrn("classif.ranger", num.trees = 500) for a classification forest.
#' @return Estimates of \eqn{Pr(D = 1 | Z)} and an 'mlr3' object of the employed model (if applicable)
#'
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

    out <- list(propensity_scores = estimator,
                mlr3_objects = NULL)

  } else if(estimator == "constant"){

    ### case 2: propensity scores are estimated by mean of D
    out <- list(propensity_scores = rep(mean(D), length(Z)),
                mlr3_objects = NULL)

  } else{

    ### case 3: propensity scores are estimated by mlr3 (or illegal input)

    # for the following choices of the learner, we require mlr3:
    if(estimator %in% c("elastic.net", "random.forest", "tree") |
       substr(estimator, start = 1, stop = 6) == "mlr3::"){

      out <- propensity_score_mlr3(Z = Z, D = D, learner = make.mlr3.environment(estimator, regr = FALSE))

    } else stop("The argument 'estimator' must be equal to either 'constant', 'elastic.net', random.forest', 'tree', an mlr3 string, or a numeric vector of the same length as Z and D!")

  } # IF

  # return
  return(out)

} # FUN



#' Estimates the baseline proxy estimator \eqn{E[Y | D=0, Z]} on the auxiliary sample
#'
#' @param Z an ( _n_ x _d_) matrix or data frame of covariates
#' @param D a binary vector of treatment status of length _n_
#' @param Y a vector of responses of length _n_
#' @param auxiliary.sample a numerical vector of indices of observations in the auxiliary sample. Length is shorter than _n_
#' @param learner the classification machine learner to be used. Either 'glm', 'random.forest', or 'tree'. Can alternatively be specified by using the mlr3 framework, for example ml_g = mlr3::lrn("regr.ranger", num.trees = 500) for a regression forest, which is also the default.
#' @param min_variation minimum variation of the predictions before random noise with distribution N(0, var(Y)/20) is added. Default is 1e-05.
#' @return Estimates of Y, both for the auxiliary sample and all observations, and an 'mlr3' object of the employed model
#'
#'
#' @export
proxy_baseline <- function(Z, D, Y,
                           auxiliary.sample,
                           learner = "random.forest",
                           min_variation = 1e-05){

  # input checks
  InputChecks_D(D)
  InputChecks_Y(Y)
  InputChecks_Z(Z)
  InputChecks_equal.length3(D, Y ,Z)

  # call main function
  proxy_baseline_NoChecks(Z = Z, D = D, Y = Y,
                          auxiliary.sample = auxiliary.sample,
                          learner = get.learner_regr(learner),
                          min_variation = min_variation)

} # END FUN


#' helper that skips the input checks
#'
#' @import mlr3 mlr3learners
#' @noRd
proxy_baseline_NoChecks <- function(Z, D, Y,
                                    auxiliary.sample,
                                    learner, # must be mlr3 object
                                    min_variation = 1e-05){

  # specify the task
  task.proxy_baseline.estimator <- mlr3::TaskRegr$new(id = "proxy_baseline",
                                                      backend = data.frame(Y, Z),
                                                      target = "Y")

  # specify that the learner predicts Y
  learner$predict_type = "response"

  # indices of the control units in the auxiliary sample
  auxiliary.sample.logical <- 1:length(Y) %in% auxiliary.sample
  idx <- which(auxiliary.sample.logical & D == 0)

  # identify the main sample
  main.sample <- setdiff(1:length(Y), auxiliary.sample)

  # fit the learner on the control units in the auxiliary sample
  learner$train(task.proxy_baseline.estimator, row_ids = idx)

  # obtain predictions for Y for all observations
  predictions.obj <- learner$predict(task.proxy_baseline.estimator)
  predictions     <- predictions.obj$response

  # if there is not much variation in the predictions, add Gaussian noise
  if(stats::var(predictions) < min_variation){

    predictions <- predictions +
      stats::rnorm(length(Y), mean = 0, sd = sqrt(stats::var(Y) / 20))

  } # IF

  # return
  return(structure(
    list(baseline.predictions.main.sample = predictions[main.sample],
         baseline.predictions.auxiliary.sample = predictions[auxiliary.sample],
         baseline.predictions.full.sample = predictions,
         auxiliary.sample = auxiliary.sample,
         mlr3_objects = list(task = task.proxy_baseline.estimator,
                             learner = learner)), class = "proxy_baseline"))

} # FUN


#' Estimates the CATE proxy estimator \eqn{E[Y | D=1, Z] - E[Y | D=0, Z]} on the auxiliary sample
#'
#' @param Z an ( _n_ x _d_) matrix or data frame of covariates
#' @param D a binary vector of treatment status of length _n_
#' @param Y a vector of responses of length _n_
#' @param auxiliary.sample a numerical vector of indices of observations in the auxiliary sample. Length is shorter than _n_
#' @param learner the regression machine learner to be used. Either 'glm', 'random.forest', or 'tree'. Can alternatively be specified by using the mlr3 framework, for example ml_g = mlr3::lrn("regr.ranger", num.trees = 500) for a regression forest, which is also the default.
#' @param proxy_baseline.estimates A vector of length _n_ of proxy estimates of the baseline estimator \eqn{E[Y | D=0, Z]}. If NULL, these will be estimated separately.
#' @param min_variation minimum variation of the predictions before random noise with distribution N(0, var(Y)/20) is added. Default is 1e-05.
#' @return Estimates of the CATE, both for the auxiliary sample and all observations, and an 'mlr3' object of each employed model
#'
#' @export
proxy_CATE <- function(Z, D, Y,
                       auxiliary.sample,
                       learner = "random.forest",
                       proxy_baseline.estimates = NULL,
                       min_variation = 1e-05){

  # input checks
  InputChecks_D(D)
  InputChecks_Y(Y)
  InputChecks_Z(Z)
  InputChecks_equal.length3(D, Y ,Z)

  # run main function
  proxy_CATE_NoChecks(Z = Z, D = D, Y = Y,
                      auxiliary.sample = auxiliary.sample,
                      learner = get.learner_regr(learner), # must be mlr3 object
                      proxy_baseline.estimates = proxy_baseline.estimates,
                      min_variation = min_variation)

} # END FUN


#' helper that skips the input checks
#'
#' @import mlr3 mlr3learners
#' @noRd
proxy_CATE_NoChecks <- function(Z, D, Y,
                                auxiliary.sample,
                                learner = "random.forest",
                                proxy_baseline.estimates = NULL,
                                min_variation = 1e-05){

  # indices of the treated units in the auxiliary sample
  auxiliary.sample.logical <- 1:length(Y) %in% auxiliary.sample
  idx.auxiliary.treated <- which(auxiliary.sample.logical & D == 1)

  # identify the main sample
  main.sample <- setdiff(1:length(Y), auxiliary.sample)

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


  if(!is.null(proxy_baseline.estimates)){

    ## if proxy baseline estimates were provided, no need for second estimation
    predictions.controls <- proxy_baseline.estimates
    mlr3.controls        <- "Not available as proxy estimates for the baseline were aprovided"

  } else{

    ## if no proxy baseline estimates were provided, estimate them
    # specify the task for estimating E[Y | D=0, Z]
    task.proxy_CATE.controls.estimator <- mlr3::TaskRegr$new(id = "cate.controls",
                                                             backend = data.frame(Y, Z),
                                                             target = "Y")

    # prepare the learner
    learner.controls <- learner

    # control units in the auxiliary sample
    idx.auxiliary.controls <- which(auxiliary.sample.logical & D == 0)

    # fit the learner on the control units in the auxiliary sample
    learner.controls$train(task.proxy_CATE.controls.estimator, row_ids = idx.auxiliary.controls)

    # obtain predictions for E[Y | D=0, Z] for all observations
    predictions.controls.obj <- learner.controls$predict(task.proxy_CATE.controls.estimator)
    predictions.controls     <- predictions.controls.obj$response

    # prepare list for output
    mlr3.controls <- list(task = task.proxy_CATE.controls.estimator,
                          learner = learner.controls)

  } # END IF

  # prepare returned objects
  predictions.Y1 <-
    list(predictions.Y1_main.sample = predictions.treated[main.sample],
         predictions.Y1_auxiliary.sample = predictions.treated[auxiliary.sample],
         predictions.Y1_full.sample = predictions.treated,
         mlr3_objects = list(task = task.proxy_CATE.treated.estimator,
                             learner = learner.treated))

  predictions.Y0 <-
    list(predictions.Y0_main.sample = predictions.controls[main.sample],
         predictions.Y0_auxiliary.sample = predictions.controls[auxiliary.sample],
         predictions.Y0_full.sample = predictions.controls,
         mlr3_objects = mlr3.controls)

  # get CATE predictions
  cate.predictions <- predictions.treated - predictions.controls

  # if there is not much variation in the predictions, add Gaussian noise
  if(stats::var(cate.predictions) < min_variation){

    cate.predictions <- cate.predictions +
      stats::rnorm(length(Y), mean = 0, sd = sqrt(stats::var(Y) / 20))

  } # IF

  return(structure(
    list(CATE.predictions.main.sample = cate.predictions[main.sample],
         CATE.predictions.auxiliary.sample = cate.predictions[auxiliary.sample],
         CATE.predictions.full.sample = cate.predictions,
         Y1.predictions = predictions.Y1,
         Y0.predictions = predictions.Y0,
         auxiliary.sample = auxiliary.sample), class = "proxy_CATE"))

} # FUN
