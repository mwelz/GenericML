#' Estimates the propensity score
#' 
#' @param Z a matrix or data frame of covariates
#' @param D a binary vector of treatment status
#' @param learner the classification machine learner to be used. Either 'glm', 'random.forest', or 'tree'. Can alternatively be specified by using the mlr3 framework, for example ml_g = mlr3::lrn("classif.ranger", num.trees = 500) for a classification forest, which is also the default. 
#' @return Estimates of Pr(D = 1 | Z) and an 'mlr3' object of the employed model 
#' 
#' @export
propensity.score <- function(Z, D, learner = "random.forest"){
  
  # input checks
  if(length(unique(D)) != 2) stop("Treatment assignment 'D' does not have 2 unique values.")
  if(!all(c(0, 1) %in% unique(D))) stop("Treatment assignment 'D' is non-binary.")
  
  # specify the task
  task.propensity.score <- mlr3::TaskClassif$new(id = "propensity.score", 
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
    
  } else{
    
    stop("Invalid argument for 'learner'. Needs to be either 'glm', 'random.forest', 'tree', or an mlr3 object")
    
  } # END IF
  
  
  # specify that the learner predicts Pr(D = 1 | Z)
  learner$predict_type = "prob"
  
  # fit the learner
  learner$train(task.propensity.score)
  
  # extract the predictions
  predictions <- learner$predict(task.propensity.score)
  
  # extract estimations of Pr(D = 1 | Z)
  probs <- predictions$prob
  
  # return 
  return(list(propensity.scores = probs[, colnames(probs) == "1"],
              mlr3.objects = list(task = task.propensity.score, 
                                  learner = learner)))
  
} # END FUN



#' Estimates the baseline proxy estimator E[Y | D=0, Z] on the auxiliary sample
#' 
#' @param Z an ( _n_ x _d_) matrix or data frame of covariates
#' @param D a binary vector of treatment status of length _n_
#' @param Y a vector of responses of length _n_
#' @param auxiliary.sample a numerical vector of indices of observations in the auxiliary sample. Length is shorter than _n_
#' @param learner the classification machine learner to be used. Either 'glm', 'random.forest', or 'tree'. Can alternatively be specified by using the mlr3 framework, for example ml_g = mlr3::lrn("regr.ranger", num.trees = 500) for a regression forest, which is also the default. 
#' @param minimum.variation minimum variation of the predictions before random noise with distribution N(0, var(Y)/20) is added. Default is 1e-05.
#' @return Estimates of Y, both for the auxiliary sample and all observations, and an 'mlr3' object of the employed model 
#' 
#' @export
baseline.proxy.estimator <- function(Z, D, Y, 
                                     auxiliary.sample,
                                     learner = "random.forest",
                                     minimum.variation = 1e-05){
  
  # input checks
  if(length(unique(D)) != 2) stop("Treatment assignment 'D' does not have 2 unique values.")
  if(!all(c(0, 1) %in% unique(D))) stop("Treatment assignment 'D' is non-binary.")
  
  # specify the task
  task.proxy.baseline.estimator <- mlr3::TaskRegr$new(id = "proxy.baseline", 
                                                      backend = data.frame(Y, Z), 
                                                      target = "Y")
  
  # specify the machine learner
  if(is.environment(learner)){
    learner <- learner
  } else if(learner == "elastic.net"){
    
    learner <- mlr3::lrn("regr.cv_glmnet", s = "lambda.min")
    
  } else if(learner == "random.forest"){
    
    learner <- mlr3::lrn("regr.ranger", num.trees = 500)
    
  } else if(learner == "tree"){
    
    learner <- mlr3::lrn("regr.rpart")
    
  } else{
    
    stop("Invalid argument for 'learner'. Needs to be either 'glm', 'random.forest', 'tree', or an mlr3 object")
    
  } # END IF
  
  # specify that the learner predicts Y
  learner$predict_type = "response"
  
  # indices of the control units in the auxiliary sample
  auxiliary.sample.logical <- 1:length(Y) %in% auxiliary.sample
  idx <- which(auxiliary.sample.logical & D == 0)
  
  # identify the main sample
  main.sample <- setdiff(1:length(Y), auxiliary.sample)
  
  # fit the learner on the control units in the auxiliary sample
  learner$train(task.proxy.baseline.estimator, row_ids = idx)
  
  # obtain predictions for Y for all observations
  predictions.obj <- learner$predict(task.proxy.baseline.estimator)
  predictions     <- predictions.obj$response
  
  # if there is not much variation in the predictions, add Gaussian noise
  if(var(predictions) < minimum.variation){
    
    predictions <- predictions +
      rnorm(length(n), mean = 0, sd = sqrt(var(Y) / 20))
    
  } # IF
  
  # return 
  return(list(baseline.predictions.main.sample = predictions[main.sample],
              baseline.predictions.auxiliary.sample = predictions[auxiliary.sample],
              baseline.predictions.full.sample = predictions,
              auxiliary.sample = auxiliary.sample,
              mlr3.objects = list(task = task.proxy.baseline.estimator, 
                                  learner = learner)))
  
} # END FUN



#' Estimates the CATE proxy estimator E[Y | D=1, Z] - E[Y | D=0, Z] on the auxiliary sample
#' 
#' @param Z an ( _n_ x _d_) matrix or data frame of covariates
#' @param D a binary vector of treatment status of length _n_
#' @param Y a vector of responses of length _n_
#' @param auxiliary.sample a numerical vector of indices of observations in the auxiliary sample. Length is shorter than _n_
#' @param learner the regression machine learner to be used. Either 'glm', 'random.forest', or 'tree'. Can alternatively be specified by using the mlr3 framework, for example ml_g = mlr3::lrn("regr.ranger", num.trees = 500) for a regression forest, which is also the default. 
#' @param proxy.baseline.estimates A vector of length _n_ of proxy estimates of the baseline estimator E[Y | D=0, Z]. If NULL, these will be estimated separately.
#' @param minimum.variation minimum variation of the predictions before random noise with distribution N(0, var(Y)/20) is added. Default is 1e-05.
#' @return Estimates of the CATE, both for the auxiliary sample and all observations, and an 'mlr3' object of each employed model 
#' 
#' @export
CATE.proxy.estimator <- function(Z, D, Y, 
                                 auxiliary.sample, 
                                 learner = "random.forest", 
                                 proxy.baseline.estimates = NULL, 
                                 minimum.variation = 1e-05){
  
  # input checks
  if(length(unique(D)) != 2) stop("Treatment assignment 'D' does not have 2 unique values.")
  if(!all(c(0, 1) %in% unique(D))) stop("Treatment assignment 'D' is non-binary.")
  
  # indices of the treated units in the auxiliary sample
  auxiliary.sample.logical <- 1:length(Y) %in% auxiliary.sample
  idx.auxiliary.treated <- which(auxiliary.sample.logical & D == 1)
  
  # identify the main sample
  main.sample <- setdiff(1:length(Y), auxiliary.sample)
  
  # specify the machine learner
  if(is.environment(learner)){
    learner <- learner
  } else if(learner == "elastic.net"){
    
    learner <- mlr3::lrn("regr.cv_glmnet", s = "lambda.min")
    
  } else if(learner == "random.forest"){
    
    learner <- mlr3::lrn("regr.ranger", num.trees = 500)
    
  } else if(learner == "tree"){
    
    learner <- mlr3::lrn("regr.rpart")
    
  } else{
    
    stop("Invalid argument for 'learner'. Needs to be either 'glm', 'random.forest', 'tree', or an mlr3 object")
    
  } # END IF
  
  # specify that the learner predicts Y
  learner$predict_type = "response"
  
  # specify the task for estimating E[Y | D=1, Z]
  task.proxy.cate.treated.estimator <- mlr3::TaskRegr$new(id = "cate.treated", 
                                                          backend = data.frame(Y, Z), 
                                                          target = "Y")
  
  # prepare the learner
  learner.treated <- learner
  
  # fit the learner on the treated units in the auxiliary sample
  learner.treated$train(task.proxy.cate.treated.estimator, row_ids = idx.auxiliary.treated)
  
  # obtain predictions for E[Y | D=1, Z] for all observations
  predictions.treated.obj <- learner.treated$predict(task.proxy.cate.treated.estimator)
  predictions.treated     <- predictions.treated.obj$response
  
  
  if(!is.null(proxy.baseline.estimates)){
    
    ## if proxy baseline estimates were provided, no need for second estimation
    predictions.controls <- proxy.baseline.estimates
    mlr3.controls        <- "Not available as proxy estimates for the baseline were aprovided"
    
  } else{
    
    ## if no proxy baseline estimates were provided, estimate them
    # specify the task for estimating E[Y | D=0, Z]
    task.proxy.cate.controls.estimator <- mlr3::TaskRegr$new(id = "cate.controls", 
                                                             backend = data.frame(Y, Z), 
                                                             target = "Y")
    
    # prepare the learner
    learner.controls <- learner
    
    # control units in the auxiliary sample
    idx.auxiliary.controls <- which(auxiliary.sample.logical & D == 0)
    
    # fit the learner on the control units in the auxiliary sample
    learner.controls$train(task.proxy.cate.controls.estimator, row_ids = idx.auxiliary.controls)
    
    # obtain predictions for E[Y | D=0, Z] for all observations
    predictions.controls.obj <- learner.controls$predict(task.proxy.cate.controls.estimator)
    predictions.controls     <- predictions.controls.obj$response
    
    # prepare list for output
    mlr3.controls <- list(task = task.proxy.cate.controls.estimator, 
                          learner = learner.controls)
    
  } # END IF
  
  # prepare returned objects
  predictions.Y1 <- 
    list(predictions.Y1_main.sample = predictions.treated[main.sample],
         predictions.Y1_auxiliary.sample = predictions.treated[auxiliary.sample],
         predictions.Y1_full.sample = predictions.treated,
         mlr3.objects = list(task = task.proxy.cate.treated.estimator,
                             learner = learner.treated))
  
  predictions.Y0 <- 
    list(predictions.Y0_main.sample = predictions.controls[main.sample],
         predictions.Y0_auxiliary.sample = predictions.controls[auxiliary.sample],
         predictions.Y0_full.sample = predictions.controls,
         mlr3.objects = mlr3.controls)
  
  # get CATE predictions
  cate.predictions <- predictions.treated - predictions.controls
  
  # if there is not much variation in the predictions, add Gaussian noise
  if(var(cate.predictions) < minimum.variation){
    
    cate.predictions <- cate.predictions +
      rnorm(length(n), mean = 0, sd = sqrt(var(Y) / 20))
    
  } # IF
  
  return(list(CATE.predictions.main.sample = cate.predictions[main.sample],
              CATE.predictions.auxiliary.sample = cate.predictions[auxiliary.sample],
              CATE.predictions.full.sample = cate.predictions,
              Y1.predictions = predictions.Y1,
              Y0.predictions = predictions.Y0,
              auxiliary.sample = auxiliary.sample))
  
} # END FUN