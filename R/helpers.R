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


#' Get the best learners across sample splits
#'
#' @param generic_targets A list of generic targets for each learner. More specifically, each component is again a list of the generic targets pertaining to the BLP, GATES, CLAN, and best analyses. Each of those lists contains a 3-dimensional array containing the generic targets of a single learner for all sample splits (except CLAN where there is one more layer of lists).
#'
#' @return A list of the names of the best learner for BLP, GATES, CLAN, and an overview which is a matrix of the lambda and lambda.bar estimates of each best learner.
#'
#' @noRd
get.best.learners <- function(generic_targets){

  learners <- names(generic_targets)

  # for each learner, take medians over the number of splits
  best.analysis <- sapply(learners,
                          function(learner) apply(generic_targets[[learner]]$best, c(1,2), stats::median))
  rownames(best.analysis) <- c("lambda", "lambda.bar")

  return(list(BLP      = learners[which.max(best.analysis["lambda", ])],
              GATES    = learners[which.max(best.analysis["lambda.bar", ])],
              CLAN     = learners[which.max(best.analysis["lambda.bar", ])],
              overview = t(best.analysis)))

} # END FUN



#' Performs VEIN
#'
#' @param generic_targets A list of generic targets for each learner. More specifically, each component is again a list of the generic targets pertaining to the BLP, GATES, CLAN, and best analyses. Each of those lists contains a 3-dimensional array containing the generic targets of a single learner for all sample splits (except CLAN where there is one more layer of lists).
#' @param best.learners.obj An object as returned by \code{get.best.learners()}, calculated on \code{generic_targets}.
#'
#' @return A list of VEIN results for all learners and only the best learners.
#'
#' @noRd
VEIN <- function(generic_targets, best.learners.obj){

  gen.ml.ls <- initialize.gen.ml(generic_targets)
  learners  <- names(generic_targets)

  for(learner in learners){

    # BLP and GATES
    for(type in c("BLP", "GATES")){

      gen.ml.ls[[type]][[learner]][,"Estimate"] <-
        apply(generic_targets[[learner]][[type]][,"Estimate",], 1, function(z) Med(z)$Med)
      gen.ml.ls[[type]][[learner]][,"CB lower"] <-
        apply(generic_targets[[learner]][[type]][,"CB lower",], 1, function(z) Med(z)$upper_median)
      gen.ml.ls[[type]][[learner]][,"CB upper"] <-
        apply(generic_targets[[learner]][[type]][,"CB upper",], 1, function(z) Med(z)$lower_median)
      p.left.raw <-
        apply(generic_targets[[learner]][[type]][,"Pr(<z)",], 1, function(z) Med(z)$lower_median)
      p.right.raw <-
        apply(generic_targets[[learner]][[type]][,"Pr(>z)",], 1, function(z) Med(z)$lower_median)
      p.left.adj  <- 2 * p.left.raw
      p.right.adj <- 2 * p.right.raw
      p.left.adj[p.left.adj > 1]   <- 1 # p values cannot exceed 1
      p.right.adj[p.right.adj > 1] <- 1
      gen.ml.ls[[type]][[learner]][,"Pr(<z) adjusted"] <- p.left.adj
      gen.ml.ls[[type]][[learner]][,"Pr(>z) adjusted"] <- p.right.adj

    } # FOR type


    # CLAN
    z.clan.nam <- names(generic_targets[[1]]$CLAN)

    for(z.clan in z.clan.nam){

      gen.ml.ls$CLAN[[learner]][[z.clan]][,"Estimate"] <-
        apply(generic_targets[[learner]][["CLAN"]][[z.clan]][,"Estimate",], 1, function(z) Med(z)$Med)
      gen.ml.ls$CLAN[[learner]][[z.clan]][,"CB lower"] <-
        apply(generic_targets[[learner]][["CLAN"]][[z.clan]][,"CB lower",], 1, function(z) Med(z)$upper_median)
      gen.ml.ls$CLAN[[learner]][[z.clan]][,"CB upper"] <-
        apply(generic_targets[[learner]][["CLAN"]][[z.clan]][,"CB upper",], 1, function(z) Med(z)$upper_median)
      p.left.raw <-
        apply(generic_targets[[learner]][["CLAN"]][[z.clan]][,"Pr(<z)",], 1, function(z) Med(z)$lower_median)
      p.right.raw <-
        apply(generic_targets[[learner]][["CLAN"]][[z.clan]][,"Pr(>z)",], 1, function(z) Med(z)$lower_median)
      p.left.adj  <- 2 * p.left.raw
      p.right.adj <- 2 * p.right.raw
      p.left.adj[p.left.adj > 1]   <- 1 # p values cannot exceed 1
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


#' helper function that returns the names of objects that shall be exported to each parallel worker via parallel::clusterExport
#' @noRd
get_varlist <- function(){
  funs <- c("get_blp.3d",
            "get_gates.3d",
            "get_best.3d",
            "get_clan.3d.ls",
            "sample_split",
            "split_fn",
            "GenericML_single_NoChecks",
            "proxy_BCA_NoChecks",
            "proxy_CATE_NoChecks",
            "setup_X1_NoChecks",
            "setup_X1_NoChecks",
            "setup_vcov_subset",
            "BLP_NoChecks",
            "BLP.classic",
            "BLP.HT",
            "get.df.from.X1_control",
            "get.vcov",
            "VEIN",
            "generic_targets_BLP",
            "quantile_group_NoChecks",
            "GATES_NoChecks", "GATES.classic", "GATES.HT", "generic_targets_GATES",
            "CLAN_NoChecks",
            "lambda_parameters_NoChecks")

  vars <- c("num.learners",
            "learners.names",
            "num.generic_targets.gates",
            "num.generic_targets.clan",
            "num.vars.in.Z_CLAN",
            "namZ_CLAN",
            "store_learners",
            "store_splits",
            "D",
            "N",
            "Z",
            "Y",
            "learners",
            "propensity_scores",
            "Z_CLAN",
            "X1_BLP",
            "X1_GATES",
            "HT",
            "vcov_BLP",
            "vcov_GATES",
            "equal_variances_CLAN",
            "quantile_cutoffs",
            "diff_GATES",
            "diff_CLAN",
            "significance_level",
            "min_variation",
            "prop_aux")
  return(c(funs, vars))
} # FUN


#' Internal function that returns the names of variables for which CLAN was performed
#' @param x A \code{"GenericML"} object
#' @noRd
CLAN_names <- function(x)
{
  return(names(x$VEIN$best_learners$CLAN))
} # FUN

# helper that checks if x is a "GenericML" object
isGenericMLcheck <- function(x)
{
  if(!inherits(x = x, what = "GenericML", which = FALSE)){
    stop("x needs to be instance of the class GenericML")
  } # IF
} # FUN


#' Internal function that returns the confidence level
#' @param x A \code{"GenericML"} object
#' @noRd
confidence_level <- function(x)
{
  return(1.0 - 2 * x$arguments$significance_level)
} # FUN
