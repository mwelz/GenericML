library(mlr3)
library(mlr3learners) # potentially a bug in mlr3; if this is not loaded, mlr::lrn won't recognize the learner class

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
#' @return Estimates of Y, both for the auxiliary sample and all observations, and an 'mlr3' object of the employed model 
#' 
#' @export
baseline.proxy.estimator <- function(Z, D, Y, 
                                     auxiliary.sample, learner = "random.forest"){
  
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
#' @return Estimates of the CATE, both for the auxiliary sample and all observations, and an 'mlr3' object of each employed model 
#' 
#' @export
CATE.proxy.estimator <- function(Z, D, Y, 
                                 auxiliary.sample, learner = "random.forest", 
                                 proxy.baseline.estimates = NULL){
  
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
  
  # return
  cate.predictions <- predictions.treated - predictions.controls
  return(list(CATE.predictions.main.sample = cate.predictions[main.sample],
              CATE.predictions.auxiliary.sample = cate.predictions[auxiliary.sample],
              CATE.predictions.full.sample = cate.predictions,
              Y1.predictions = predictions.Y1,
              Y0.predictions = predictions.Y0,
              auxiliary.sample = auxiliary.sample))
  
} # END FUN


#' Performs generic ML for a given learning technique (with only one split of the data)
#' 
#' @param Z a matrix or data frame of covariates
#' @param D a binary vector of treatment status of length 
#' @param Y a vector of responses of length
#' @param propensity.scores a vector of propensity scores
#' @param learner The machine learner that shall be used
#' @param M.set main set
#' @param A.set auxiliary set
#' @param quantile.cutoffs Cutoff points of quantiles that shall be used for GATES grouping
#' 
#' TODO: instructions on how mlr3 input is supposed to work (needs to be a string!)
#' TODO: comments on CLAN: If there are categorical variables, apply one-hot-encoding to Z.clan. The interpretation then becomes: Is there a factor that is overproportionally present in the least or most affected group?
#' 
#' @export
get.generic.ml.for.given.learner <- function(Z, D, Y, 
                                             propensity.scores,
                                             learner = 'mlr3::lrn("cv_glmnet", s = "lambda.min")',
                                             M.set, A.set,
                                             Z.clan = NULL, 
                                             proportion.in.main.set = 0.5, 
                                             quantile.cutoffs = c(0.25, 0.5, 0.75),
                                             significance.level = 0.05){
  
  ### step 1: input checks ---- 
  if(is.null(Z.clan)) Z.clan <- Z # if no input provided, set it equal to Z
  
  ### step 2a: learn proxy predictors by using the auxiliary set ----
  
  # get the proxy baseline estimator for the main sample
  proxy.baseline.obj <- baseline.proxy.estimator(Z = Z, D = D, Y = Y, 
                                                 auxiliary.sample = A.set, 
                                                 learner = make.mlr3.string(learner, regr = TRUE))
  proxy.baseline     <- proxy.baseline.obj$baseline.predictions.main.sample
  
  # get the proxy estimator of the CATE for the main sample
  proxy.cate.obj <- 
    CATE.proxy.estimator(Z = Z, D = D, Y = Y,
                         auxiliary.sample = A.set, 
                         learner = make.mlr3.string(learner, regr = TRUE),
                         proxy.baseline.estimates = proxy.baseline.obj$baseline.predictions.full.sample)
  proxy.cate <- proxy.cate.obj$CATE.predictions.main.sample
  
  
  ### step 2b: estimate BLP parameters by OLS (TODO: HT transformation!) ----
  blp.obj <- get.BLP.params.classic(D = D[M.set], Y = Y[M.set],
                                    propensity.scores = propensity.scores[M.set],
                                    proxy.baseline = proxy.baseline, 
                                    proxy.cate = proxy.cate, 
                                    significance.level = significance.level)
  
  
  ### step 2c: estimate GATES parameters by OLS (TODO: HT transformation!) ----
  # group the proxy estimators for the CATE in the main sample by quantiles. TODO: intervals need to be [) instead of (]
  group.membership.main.sample <- quantile.group(proxy.cate, 
                                                 cutoffs = quantile.cutoffs, 
                                                 quantile.nam = TRUE) 
  
  gates.obj <- get.GATES.params.classic(D = D[M.set], Y = Y[M.set],
                                        propensity.scores = propensity.scores[M.set],
                                        group.membership.main.sample = group.membership.main.sample, 
                                        proxy.baseline = proxy.baseline, proxy.cate = proxy.cate,
                                        significance.level = significance.level)
  
  
  ### step 2d: estimate CLAN parameters in the main sample
  clan.obj <- get.CLAN.parameters(Z.clan.main.sample = Z.clan[M.set,], 
                                  group.membership.main.sample = group.membership.main.sample)
  
  
  ### step 2e: get parameters over which we maximize to find the "best" ML method ----
  best.obj <- best.ml.method.parameters(BLP.obj = blp.obj, 
                                        GATES.obj = gates.obj, 
                                        proxy.cate.main.sample = proxy.cate,
                                        group.membership.main.sample = group.membership.main.sample)
  
  ### organize output in a list ----
  return(list(BLP = blp.obj, 
              GATES = gates.obj, 
              CLAN = clan.obj,
              best = best.obj,
              CATE.proxy = proxy.cate.obj,
              baseline.proxy = proxy.baseline,
              group.membership_M.set = group.membership.main.sample))
  
} # END FUN



VEIN <- function(generic.ml.across.learners.obj, best.learners.obj){
  
  gen.ml.ls <- initialize.gen.ml(generic.ml.across.learners.obj)
  learners  <- names(generic.ml.across.learners.obj)
  
  for(learner in learners){
    
    # BLP and GATES
    for(type in c("BLP", "GATES")){
      
      gen.ml.ls[[type]][[learner]][,"Estimate"] <- 
        apply(generic.ml.across.learners.obj[[learner]][[type]][,"Estimate",], 1, function(z) Med(z)$Med)
      gen.ml.ls[[type]][[learner]][,"CB lower"] <- 
        apply(generic.ml.across.learners.obj[[learner]][[type]][,"CB lower",], 1, function(z) Med(z)$upper.median)
      gen.ml.ls[[type]][[learner]][,"CB upper"] <- 
        apply(generic.ml.across.learners.obj[[learner]][[type]][,"CB upper",], 1, function(z) Med(z)$lower.median)
      pval <- 
        apply(generic.ml.across.learners.obj[[learner]][[type]][,"Pr(>|z|)",], 1, function(z) Med(z)$lower.median)
      gen.ml.ls[[type]][[learner]][,"p-value raw"] <- pval
      pval.adj <- 2 * pval
      pval.adj[pval.adj > 1] <- 1 # p-values cannot exceed 1
      gen.ml.ls[[type]][[learner]][,"p-value adjusted"] <- pval.adj
      
    } # FOR type
    
    
    # CLAN
    z.clan.nam <- names(generic.ml.across.learners.obj[[1]]$CLAN)
    
    for(z.clan in z.clan.nam){
      
      gen.ml.ls$CLAN[[learner]][[z.clan]][,"Estimate"] <- 
        apply(generic.ml.across.learners.obj[[learner]][["CLAN"]][[z.clan]][,"Estimate",], 1, function(z) Med(z)$Med)
      gen.ml.ls$CLAN[[learner]][[z.clan]][,"CB lower"] <- 
        apply(generic.ml.across.learners.obj[[learner]][["CLAN"]][[z.clan]][,"CB lower",], 1, function(z) Med(z)$upper.median)
      gen.ml.ls$CLAN[[learner]][[z.clan]][,"CB upper"] <- 
        apply(generic.ml.across.learners.obj[[learner]][["CLAN"]][[z.clan]][,"CB upper",], 1, function(z) Med(z)$upper.median)
      pval <- 
        apply(generic.ml.across.learners.obj[[learner]][["CLAN"]][[z.clan]][,"Pr(>|z|)",], 1, function(z) Med(z)$lower.median)
      gen.ml.ls$CLAN[[learner]][[z.clan]][,"p-value raw"] <- pval
      pval.adj <- 2 * pval
      pval.adj[pval.adj > 1] <- 1 # p-values cannot exceed 1
      gen.ml.ls$CLAN[[learner]][[z.clan]][,"p-value adjusted"] <- pval.adj
      
    } # FOR Z.clan
    
  } # FOR learners
  

  return(list(best.learners = list(BLP = gen.ml.ls$BLP[[best.learners.obj$best.learner.for.CATE]],
                                   GATES = gen.ml.ls$GATES[[best.learners.obj$best.learner.for.GATES]],
                                   CLAN = gen.ml.ls$CLAN[[best.learners.obj$best.learner.for.GATES]]),
              all.learners = gen.ml.ls))
  
} # END FUN


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
#' @param proportion.in.main.set proportion of samples that shall be in main set. Default is 0.5.
#' @param significance.level significance level for VEIN. Default is 0.05.
#' @param store.learners Logical. If TRUE, all intermediate results of the learners will be stored. Default is FALSE.
#' @param store.splits Logical. If `TRUE`, information on the sample splits will be stored. Default is `FALSE`.
#' 
#' @export
genericML <- function(Z, D, Y, 
                      learner.propensity.score = "elastic.net", 
                      learners.genericML,
                      num.splits = 100,
                      Z.clan = NULL,
                      quantile.cutoffs = c(0.25, 0.5, 0.75),
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
