make.mlr3.environment <- function(learner.str, regr = TRUE){
  # helper function. Requires input of type 'mlr3::lrn("cv_glmnet", s = "lambda.min")' (note the absence of classif and regr)

  if(substr(learner.str, start = 1, stop = 6) != "mlr3::"){

    learner <- learner.str

  } else{

    learner <- paste0(substr(learner.str, start = 1, stop = 11), ifelse(regr, "regr.", "classif."),
                      substr(learner.str, start = 12, stop = 1e8))

    learner <- eval(parse(text = learner))

  } # END IF

  return(learner)

} # END FUN



get.best.learners <- function(generic.ml.across.learners.obj){

  learners <- names(generic.ml.across.learners.obj)

  # for each learner, take medians over the number of splits
  best.analysis <- sapply(learners,
                          function(learner) apply(generic.ml.across.learners.obj[[learner]]$best, c(1,2), stats::median))
  rownames(best.analysis) <- c("lambda", "lambda.bar")

  return(list(BLP      = learners[which.max(best.analysis["lambda", ])],
              GATES    = learners[which.max(best.analysis["lambda.bar", ])],
              CLAN     = learners[which.max(best.analysis["lambda.bar", ])],
              overview = t(best.analysis)))

} # END FUN



# TODO: make documentation; generic.ml.across.learners.obj is NOT the correct description for the input anymore (check genericML function!)
VEIN <- function(generic.ml.across.learners.obj, best.learners.obj){

  gen.ml.ls <- initialize.gen.ml(generic.ml.across.learners.obj)
  learners  <- names(generic.ml.across.learners.obj)

  for(learner in learners){

    # BLP and GATES
    for(type in c("BLP", "GATES")){

      gen.ml.ls[[type]][[learner]][,"Estimate"] <-
        apply(generic.ml.across.learners.obj[[learner]][[type]][,"Estimate",], 1, function(z) Med(z)$Med)
      gen.ml.ls[[type]][[learner]][,"CB lower"] <-
        apply(generic.ml.across.learners.obj[[learner]][[type]][,"CB lower",], 1, function(z) Med(z)$upper_median)
      gen.ml.ls[[type]][[learner]][,"CB upper"] <-
        apply(generic.ml.across.learners.obj[[learner]][[type]][,"CB upper",], 1, function(z) Med(z)$lower_median)
      p.left.raw <-
        apply(generic.ml.across.learners.obj[[learner]][[type]][,"Pr(<z)",], 1, function(z) Med(z)$lower_median)
      p.right.raw <-
        apply(generic.ml.across.learners.obj[[learner]][[type]][,"Pr(>z)",], 1, function(z) Med(z)$lower_median)
      p.left.adj  <- 2 * p.left.raw
      p.right.adj <- 2 * p.right.raw
      p.left.adj[p.left.adj > 1]   <- 1 # p-values cannot exceed 1
      p.right.adj[p.right.adj > 1] <- 1
      gen.ml.ls[[type]][[learner]][,"Pr(<z) adjusted"] <- p.left.adj
      gen.ml.ls[[type]][[learner]][,"Pr(>z) adjusted"] <- p.right.adj

    } # FOR type


    # CLAN
    z.clan.nam <- names(generic.ml.across.learners.obj[[1]]$CLAN)

    for(z.clan in z.clan.nam){

      gen.ml.ls$CLAN[[learner]][[z.clan]][,"Estimate"] <-
        apply(generic.ml.across.learners.obj[[learner]][["CLAN"]][[z.clan]][,"Estimate",], 1, function(z) Med(z)$Med)
      gen.ml.ls$CLAN[[learner]][[z.clan]][,"CB lower"] <-
        apply(generic.ml.across.learners.obj[[learner]][["CLAN"]][[z.clan]][,"CB lower",], 1, function(z) Med(z)$upper_median)
      gen.ml.ls$CLAN[[learner]][[z.clan]][,"CB upper"] <-
        apply(generic.ml.across.learners.obj[[learner]][["CLAN"]][[z.clan]][,"CB upper",], 1, function(z) Med(z)$upper_median)
      p.left.raw <-
        apply(generic.ml.across.learners.obj[[learner]][["CLAN"]][[z.clan]][,"Pr(<z)",], 1, function(z) Med(z)$lower_median)
      p.right.raw <-
        apply(generic.ml.across.learners.obj[[learner]][["CLAN"]][[z.clan]][,"Pr(>z)",], 1, function(z) Med(z)$lower_median)
      p.left.adj  <- 2 * p.left.raw
      p.right.adj <- 2 * p.right.raw
      p.left.adj[p.left.adj > 1]   <- 1 # p-values cannot exceed 1
      p.right.adj[p.right.adj > 1] <- 1
      gen.ml.ls$CLAN[[learner]][[z.clan]][,"Pr(<z) adjusted"] <- p.left.adj
      gen.ml.ls$CLAN[[learner]][[z.clan]][,"Pr(>z) adjusted"] <- p.right.adj

    } # FOR Z_CLAN

  } # FOR learners


  return(list(best_learners = list(BLP = gen.ml.ls$BLP[[best.learners.obj$BLP]],
                                   GATES = gen.ml.ls$GATES[[best.learners.obj$GATES]],
                                   CLAN = gen.ml.ls$CLAN[[best.learners.obj$CLAN]]),
              all_learners = gen.ml.ls))

} # END FUN


#' helper function that calculates an error covariance matrix estimator of a linear model
#'
#' (internal use only)
#'
#' @param x a linear model object
#' @param vcov_control a list with two elements called 'estimator' and 'arguments'. The argument 'estimator' is a string specifying the covariance matrix estimator to be used; specifies a covariance estimator function in the sandwich package (https://cran.r-project.org/web/packages/sandwich/sandwich.pdf). Recommended estimators are "vcovBS", "vcovCL", "vcovHAC", and "vcovHC". Default is 'vcovHC'. The element 'arguments' is a list of arguments that shall be passed to the function specified in the element 'estimator'. Default leads to the (homoskedastic) ordinary least squares covariance matrix estimator. See the reference manual of the sandwich package for details (https://cran.r-project.org/web/packages/sandwich/vignettes/sandwich.pdf).
#'
#' @import sandwich
#'
#' @noRd
get.vcov <- function(x,
                     vcov_control = setup_vcov()){

  # append the model so that we can pass this list to do.call
  arguments   <- vcov_control$arguments
  arguments$x <- x

  # return the estimate
  do.call(what = eval(parse(text = paste0("sandwich::", vcov_control$estimator))),
          args = arguments)
  # do.call(what = vcov_control$estimator, args = arguments)

} # FUN


# helper function to prepare the custom part of the regressor matrix in BLP and GATES
get.df.from.X1_control <- function(functions.of.Z_mat,
                                     X1_control){

  custom        <- X1_control$covariates
  fixed.eff     <- X1_control$fixed_effects
  out           <- data.frame(functions.of.Z_mat[, X1_control$funs_Z, drop = FALSE])
  colnames(out) <- X1_control$funs_Z

  if(!is.null(fixed.eff)) out$fixed.effects <- factor(fixed.eff)
  if(!is.null(custom))    out <- data.frame(out, custom)

  out

} # FUN


#' Performs the sample splitting (internal use)
#'
#' @param D Binary vector of treatment assignment
#' @param N sample size
#' @param N_set 1:N
#' @param prop Total number of samples in the auxiliary set
#'
#' @noRd
sample_split <- function(D, N, N_set = 1:N, prop){

  temp <- TRUE

  while(temp){

    # sample candidate set for A_set
    A_set <- sample(N_set, size = prop, replace = FALSE)

    # Avoid imbalance in A_set for estimation of BCA and CATE.
    # BCA is estimated on the control units, CATE on the treated units.
    # We want to avoid that either of them is estimated on too small a sample.
    # To achieve this, have the control units in A_set make up no more
    # than 90% of all samples in A_set
    if(mean(D[A_set] == 0) <= 0.9) temp <- FALSE

  } # WHILE

  A_set <- as.integer(sort(A_set, decreasing = FALSE))

  # return
  return(list(A_set = A_set,
              M_set = setdiff(N_set, A_set)))

} # FUN
