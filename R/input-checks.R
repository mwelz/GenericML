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
  if((!all(c(0, 1) %in% unique(D))) | (length(unique(D)) != 2)) stop("Treatment assignment D is not binary", call. = FALSE)
  if(!is.numeric(D)) stop("D is not numeric", call. = FALSE)

} # FUN


InputChecks_Y <- function(Y){

  if(!is.numeric(Y)) stop("Y is not numeric", call. = FALSE)
  if(!(is.matrix(Y) | is.data.frame(Y) | is.vector(Y))) stop("Y must be either a vector, data frame, or matrix", call. = FALSE)

  if(is.matrix(Y) | is.data.frame(Y)){

    if(ncol(Y) != 1) stop("The matrix/data frame Y must have only one column", call. = FALSE)

  } # IF

} # FUN


InputChecks_Z <- function(Z){

  if(!is.numeric(Z)) stop("Z is not numeric", call. = FALSE)
  if(!(is.matrix(Z) | is.data.frame(Z))) stop("Z must be either a data frame or matrix", call. = FALSE)

} # FUN


InputChecks_Z_CLAN <- function(Z_CLAN){

  if(!is.null(Z_CLAN)){

    if(!(is.matrix(Z_CLAN) | is.data.frame(Z_CLAN) )){
      stop("Z_CLAN must be either a data frame or matrix", call. = FALSE)
    }
  }

} # FUN




# helper that throws error in case of illegal input in 'X1_control'
InputChecks_X1 <- function(X1_control, num.obs){

  if(!all(c("functions_of_Z", "custom_covariates", "fixed_effects") %in% names(X1_control))){
      stop(paste0("The list ", deparse(substitute(X1_control)), "must consist of three elements called ",
                  "'functions_of_Z', 'custom_covariates', and 'fixed_effects"), call. = FALSE)
    } # IF


  if(!(is.matrix(X1_control$custom_covariates) | is.vector(X1_control$custom_covariates) | is.null(X1_control$custom_covariates))){
    stop(paste0(deparse(substitute(X1_control)), "$custom_covariates must be either NULL, a vector, or a matrix"), call. = FALSE)
  } # IF


  if(!is.null(X1_control$custom_covariates)){
    if(!is.numeric(X1_control$custom_covariates)) stop(paste0(deparse(substitute(X1_control)), "$custom_covariates must be numeric or NULL"), call. = FALSE)

    if(is.vector(X1_control$custom_covariates)){

      if(length(X1_control$custom_covariates) != num.obs){
        stop(paste0(deparse(substitute(X1_control)), "$custom_covariates must be NULL or of the same length as Y"), call. = FALSE)
      } # IF

    } else{

      if(nrow(X1_control$custom_covariates) != num.obs){
        stop(paste0(deparse(substitute(X1_control)), "$custom_covariates must be NULL or of the same length as Y"), call. = FALSE)
      } # IF

    } # IF

  } # IF


  if(!(is.vector(X1_control$fixed_effects) | is.null(X1_control$fixed_effects))){
    stop(paste0(deparse(substitute(X1_control)), "$fixed_effects must be either NULL or a vector"), call. = FALSE)
  } # IF


  if(is.vector(X1_control$fixed_effects)){

    if(!is.numeric(X1_control$fixed_effects)){
      stop(paste0(deparse(substitute(X1_control)), "$fixed_effects must be numeric"), call. = FALSE)
    } else{

      if(length(X1_control$fixed_effects) != num.obs){
        stop(paste0(deparse(substitute(X1_control)), "$fixed_effects must be NULL or of the same length as Y"),
             call. = FALSE)
      } # IF

    }# IF
  } # IF


  legalinput <- X1_control$functions_of_Z %in% c("S", "B", "p")

  if(!all(legalinput)){

    stop(paste0("Entries '",
                paste(X1_control$functions_of_Z[!legalinput], collapse = "', '"),
                "' of ", deparse(substitute(X1_control)),
                "$functions_of_Z are not contained in c('S', 'B', 'p')!"), call. = FALSE)

  } # IF

  if(!is.null(X1_control$fixed_effects) & !is.vector(X1_control$fixed_effects)){
    stop("The fixed effects need to be a vector", call. = FALSE)
  } # IF

} # FUN



InputChecks_vcov.control <- function(vcov_control){

  if(!is.list(vcov_control)) stop(paste0(deparse(substitute(vcov_control))),
                                  " must be a list", call. = FALSE)

  if(length(vcov_control) != 2) stop(paste0("The list ", deparse(substitute(vcov_control))),
                                     " must be of length 2", call. = FALSE)

  if(!all(c("estimator", "arguments") %in% names(vcov_control))){

    stop(paste0("The list ", deparse(substitute(vcov_control))),
         " must have two elements called 'estimator' and 'arguments', respectively", call. = FALSE)

  } # IF


  if(!vcov_control$estimator %in% c("vcovBS", "vcovCL", "vcovHAC", "vcovHC")){
    stop(paste0("The element ", deparse(substitute(vcov_control))), "$estimator",
         " needs to be in c('vcovBS', 'vcovCL', 'vcovHAC', 'vcovHC')", call. = FALSE)
  } # IF

  if(!is.list(vcov_control$arguments)){

    stop(paste0(deparse(substitute(vcov_control))),
         "$arguments must be a list", call. = FALSE)

  } else{

    if(!"type" %in% names(vcov_control$arguments)) stop(paste0("The list ", deparse(substitute(vcov_control))),
                                                 "$arguments must contain an element called 'type'",
                                                 call. = FALSE)

  } # IF

} # FUN


InputChecks_diff <- function(diff, K){

  if(!is.list(diff)) stop(paste0(deparse(substitute(diff))),
                                         " must be a list", call. = FALSE)

  if(length(diff) != 2) stop(paste0("The list ", deparse(substitute(diff))),
                                       " must be of length 2", call. = FALSE)

  if(!all(c("subtract_from", "subtracted") %in% names(diff))){

    stop(paste0("The list ", deparse(substitute(diff))),
         " must have two elements called 'subtract_from' and 'subtracted', respectively", call. = FALSE)

  } # IF


  if(!diff$subtract_from %in% c("most", "least")){
    stop(paste0("The element ", deparse(substitute(diff)), "$subtract_from",
         " must be equal to 'most' or 'least'"), call. = FALSE)
  }

  if(!(is.vector(diff$subtracted) | is.numeric(diff$subtracted))){

    stop(paste0(deparse(substitute(diff)), "$subtracted",
                " must be a numeric vector"), call. = FALSE)

  }

  if(any(diff$subtracted < 1) | any(diff$subtracted > K)){
    stop(paste0("The numeric vector ", deparse(substitute(diff)), "$subtracted",
                " must be a subset of {1,2,...,K}, where K = ", K, " is the number of groups with the supplied arguments (controlled through the argument 'quantile_cutoffs')"), call. = FALSE)
  }

} # FUN


InputChecks_group.membership <- function(group.membership){

  if(is.null(attr(group.membership, which = "type"))) stop(paste0("The object ",
                                                       deparse(substitute(group.membership)),
                                                       " needs to be returned by quantile_group()"))

  if(attr(group.membership, which = "type") != "quantile_group") stop(paste0("The object ",
                                                                  deparse(substitute(group.membership)),
                                                                  " needs to be returned by quantile_group()"))


} # FUN


InputChecks_proxy.estimators <- function(proxy.estimators, baseline = TRUE){

  if(baseline){

    if(!inherits(proxy.estimators, what = "proxy_baseline")){

      stop(paste0("The object ",
                  deparse(substitute(proxy.estimators)),
                  " needs to be an instance of proxy_baseline()"))
    }

  } else{

    if(!inherits(proxy.estimators, what = "proxy_CATE")){

      stop(paste0("The object ",
                  deparse(substitute(proxy.estimators)),
                  " needs to be an instance of proxy_CATE()"))
    }

  }

} # FUN


# checks if learner is correctly specified. If yes, that learner is returned
get.learner_regr <- function(learner){

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

    stop("Invalid argument for 'learner'. Needs to be either 'elastic.net', 'random.forest', 'tree', or an mlr3 object")

  } # END IF

  return(learner)

} # FUN
