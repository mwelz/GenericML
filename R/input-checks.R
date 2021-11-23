# void functions for input checks

InputChecks_equal.length3 <- function(x, y, z){

  X <- as.matrix(x)
  Y <- as.matrix(y)
  Z <- as.matrix(z)

  if(!(nrow(X) == nrow(Y) & nrow(X) == nrow(Z) & nrow(Z) == nrow(Y))){
    stop(paste0(deparse(substitute(x)), ", ",
                deparse(substitute(y)), ", ",
                deparse(substitute(z)),
                " need to have an equal number of observations"), call. = FALSE)
  }
} # FUN

InputChecks_equal.length2 <- function(x, y){

  X <- as.matrix(x)
  Y <- as.matrix(y)

  if(!(nrow(X) == nrow(Y))){
    stop(paste0(deparse(substitute(x)), ", ",
                deparse(substitute(y)),
                " need to have an equal number of observations"), call. = FALSE)
  }
} # FUN



InputChecks_D <- function(D){

  # input checks
  if(any(is.na(D))) stop("D contains missing values", call. = FALSE)
  if(!(is.numeric(D) & is.vector(D))) stop("D must be a numeric vector", call. = FALSE)
  if((!all(c(0, 1) %in% unique(D))) | (length(unique(D)) != 2)) stop("Treatment assignment D is not binary", call. = FALSE)

} # FUN


InputChecks_Y <- function(Y){

  # input checks
  if(any(is.na(Y))) stop("Y contains missing values", call. = FALSE)
  if(!(is.numeric(Y) & is.vector(Y))) stop("Y must be a numeric vector", call. = FALSE)

} # FUN


InputChecks_Z <- function(Z){

  # input checks
  if(any(is.na(Z))) stop("Z contains missing values", call. = FALSE)
  if(!(is.numeric(Z) & is.matrix(Z))) stop("Z must be a numeric matrix. Did you supply a data frame?", call. = FALSE)

} # FUN


InputChecks_Z_CLAN <- function(Z_CLAN){

  if(!is.null(Z_CLAN)){

    # input checks
    if(any(is.na(Z_CLAN))) stop("Z_CLAN contains missing values", call. = FALSE)
    if(!(is.numeric(Z_CLAN) & is.matrix(Z_CLAN))) stop("Z_CLAN must be a numeric matrix or NULL. Did you supply a data frame?", call. = FALSE)

  } # IF

} # FUN




# helper that throws error in case of illegal input in 'X1_control'
InputChecks_X1 <- function(X1_control, num.obs){

  if(class(X1_control) != "setup_X1"){
    stop(paste0(deparse(substitute(X1_control))),
         " must be an instance of setup_X1()", call. = FALSE)
  } # IF


  if(!is.null(X1_control$covariates)){

    if(nrow(X1_control$covariates) != num.obs){
      stop(paste0(deparse(substitute(X1_control)),
                  "$covariates must have the same number of rows as Z"),
           call. = FALSE)
    } # IF
  } # IF !NULL


  if(!is.null(X1_control$fixed_effects)){

    if(length(X1_control$fixed_effects) != num.obs){
      stop(paste0(deparse(substitute(X1_control)),
                  "$fixed_effects must have the same length as Y"),
           call. = FALSE)
    } # IF
  } # IF !NULL

} # FUN



InputChecks_vcov.control <- function(vcov_control){

  if(class(vcov_control) != "setup_vcov"){
    stop(paste0(deparse(substitute(vcov_control))),
                " must be an instance of setup_vcov()", call. = FALSE)
  } # IF

} # FUN


InputChecks_diff <- function(diff, K){

  if(class(diff) != "setup_diff"){
    stop(paste0(deparse(substitute(diff))),
         " must be an instance of setup_diff()", call. = FALSE)
  } # IF

  if(any(diff$subtracted < 1) | any(diff$subtracted > K)){
    stop(paste0("The numeric vector ", deparse(substitute(diff)), "$subtracted",
                " must be a subset of {1,2,...,K}, where K = ", K, " is the number of groups that were supplied  (controlled through the argument 'quantile_cutoffs'). K is equal to the cardinality of 'quantile_cutoffs' plus one."), call. = FALSE)
  }

  if(diff$subtract_from == "most" & K %in% diff$subtracted){
    stop("The most affected group cannot be subtracted from itself")
  }

  if(diff$subtract_from == "least" & 1 %in% diff$subtracted){
    stop("The least affected group cannot be subtracted from itself")
  }

} # FUN


InputChecks_group.membership <- function(group.membership){

  if(is.null(attr(group.membership, which = "type"))) stop(paste0("The object ",
                                                                  deparse(substitute(group.membership)),
                                                                  " must be returned by quantile_group()"))

  if(attr(group.membership, which = "type") != "quantile_group") stop(paste0("The object ",
                                                                             deparse(substitute(group.membership)),
                                                                             " must be returned by quantile_group()"))


} # FUN


InputChecks_proxy.estimators <- function(proxy.estimators, baseline = TRUE){

  if(baseline){

    if(!inherits(proxy.estimators, what = "proxy_BCA")){

      stop(paste0("The object ",
                  deparse(substitute(proxy.estimators)),
                  " needs to be an instance of proxy_BCA()"))
    }

  } else{

    if(!inherits(proxy.estimators, what = "proxy_CATE")){

      stop(paste0("The object ",
                  deparse(substitute(proxy.estimators)),
                  " needs to be an instance of proxy_CATE()"))
    }

  }

} # FUN


InputChecks_num_splits <- function(num_splits){
  stopifnot(length(num_splits) == 1)
  stopifnot(num_splits %% 1 == 0)
  if(num_splits < 2){
    stop(paste0("num_splits must be equal to at least 2. If you want to run GenericML() for ",
                "a single split, please use the function GenericML_single()."), call. = FALSE)
  }
} # FUN



# checks if learner is correctly specified. If yes, that learner is returned
get.learner_regr <- function(learner){

  # specify the machine learner
  if(is.environment(learner)){
    learner <- learner
  } else if(learner == "lasso"){

    learner <- mlr3::lrn("regr.cv_glmnet", s = "lambda.min", alpha = 1)

  } else if(learner == "random_forest"){

    learner <- mlr3::lrn("regr.ranger", num.trees = 500)

  } else if(learner == "tree"){

    learner <- mlr3::lrn("regr.rpart")

  } else{

    stop("Invalid argument for 'learner'. Needs to be either 'lasso', 'random_forest', 'tree', or an mlr3 object")

  } # END IF

  return(learner)

} # FUN


InputChecks_propensity_scores <- function(propensity_scores){

  # check if data is from a randomized experiment
  if(any(propensity_scores > 0.65 | propensity_scores < 0.35)){
    warning(paste0("Some propensity scores are outside the ",
                   "interval [0.35, 0.65]. In a randomized experiment, we would ",
                   "expect all propensity scores to be equal to roughly 0.5. ",
                   "The theory of the paper ",
                   "is only valid for randomized experiments. Are ",
                   "you sure your data is from a randomomized experiment ",
                   "and the estimator of the scores has been chosen appropriately?"),
            call. = FALSE)
  } # IF

  if(any(propensity_scores > 0.95)){

    stop(paste0("Some estimated propensity scores are higher than 0.95, ",
                " which is not sufficiently bounded away from one.",
                " Are you sure your data is from a randomomized experiment ",
                "and the estimator of the scores has been chosen appropriately?"))
  } # IF

  if(any(propensity_scores < 0.05)){

    stop(paste0("Some estimated propensity scores are lower than 0.05, ",
                " which is not sufficiently bounded away from zero.",
                " Are you sure your data is from a randomomized experiment ",
                "and the estimator of the scores has been chosen appropriately?"))
  } # IF

} # FUN


InputChecks_index_set <- function(set, num_obs){

  stopifnot(is.numeric(set) & is.vector(set))

  if(any(set %% 1 != 0)){
    stop("The indices in the index set must be index-valued")
  }

  if(min(set) < 0){
    stop("The indices in the index set must be strictly positive")
  }

  if(max(set) > num_obs){
    stop("The largest index in the index set cannot be larger than the number of observations")
  }

  if(any(duplicated(set))){
    stop("All indices in the index set must be unique")
  }

} # FUN


#' Check if user's OS is a Unix system
#'
#' @return
#' A Boolean that is \code{TRUE} if the user's operating system is a Unix system and \code{FALSE} otherwise.
#'
#' @export
TrueIfUnix <- function(){
  .Platform$OS.type == "unix"
}
