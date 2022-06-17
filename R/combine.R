#' Combine several GenericML objects
#'
#' This function combines multiple \code{"\link{GenericML}"} objects into one  \code{"\link{GenericML}"} object. Combining several  \code{"\link{GenericML}"} objects can be useful when you cannot run \code{\link{GenericML}()} for sufficiently many splits due to memory constraints. In this case, you may run \code{\link{GenericML}()} multiple times with only a small number of sample splits each and combine the returned \code{"\link{GenericML}"} objects into one \code{GenericML} object with this function.
#'
#' @param x A list of \code{"\link{GenericML}"} objects, as returned by the function \code{\link{GenericML}()}.
#'
#' @return
#' A\code{"\link{GenericML}"} object as returned by \code{\link{GenericML}()}. In the \code{arguments} component of this object, the objects \code{parallel}, \code{num_cores}, \code{seed}, and \code{store_learners} are set to \code{NULL} as these might differ between the individual \code{GenericML} objects in \code{x}. Moreover, the \code{propensity_scores} component of the returned object is taken from the first \code{"\link{GenericML}"} object in \code{x}.
#'
#' @details
#' To ensure consistency of the estimates, all \code{"\link{GenericML}"} objects in the list \code{x} must have the exact same parameter specifications in their original call to \code{\link{GenericML}()}, except for the parameters \code{num_splits}, \code{parallel}, \code{num_cores}, \code{seed}, and \code{store_learners} (i.e. these arguments may vary between the \code{"\link{GenericML}"} objects in the list \code{x}). An error will be thrown if this is not satisfied.
#'
#' @seealso
#' \code{\link{GenericML}()}
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
#' ## call GenericML three times and store the returned objects in a list x
#' x <- lapply(1:3, function(...) GenericML(Z, D, Y,
#'                                learners, num_splits = 2,
#'                                parallel = FALSE))
#'
#' ## combine the objects in x into one GenericML object
#' genML <- GenericML_combine(x)
#'
#' ## you can use all methods of GenericML objects on the combined object, for instance accessors:
#' get_BLP(genML, plot = TRUE)
#' }
#'
#' @export
GenericML_combine <- function(x)
{
  ## input checks
  InputChecks_GenericML_combine(x)

  ## function that binds 3d arrays along the third dimension
  f <-  function(y, z) abind::abind(y, z, along = 3L)

  ## initialize
  m               <- length(x)
  vars            <- names(x[[1]]$generic_targets[[1]]$CLAN)
  learners        <- names(x[[1]]$generic_targets)
  generic_targets <- rep(list(NULL), length(learners))


  ## merge the generic target 3D arrays
  for(i in seq_along(learners))
  {

    blp_ls   <- lapply(1:m, function(j) x[[j]]$generic_targets[[i]]$BLP)
    gates_ls <- lapply(1:m, function(j) x[[j]]$generic_targets[[i]]$GATES)
    best_ls  <- lapply(1:m, function(j) x[[j]]$generic_targets[[i]]$best)

    clan <- list()
    for(k in seq_along(vars)){
      clan_ls  <- lapply(1:m, function(j) x[[j]]$generic_targets[[i]]$CLAN[[k]])
      clan[[k]] <- Reduce(f = f, x = clan_ls)
    } # FOR
    names(clan) <- vars

    generic_targets[[i]] <- list(BLP = Reduce(f = f, x = blp_ls),
                                 GATES = Reduce(f = f, x = gates_ls),
                                 CLAN = clan,
                                 best = Reduce(f = f, x = best_ls))

  } # FOR i
  names(generic_targets) <- learners

  ## find the best learners
  best.learners <- get.best.learners(generic_targets = generic_targets)

  ## perform VEIN analysis
  vein <- VEIN(generic_targets = generic_targets, best.learners.obj = best.learners)

  # prepare split output
  if(!is.null(x[[1]]$arguments$store_splits))
  {
    splits <- Reduce(f = cbind, x = lapply(1:m, function(i) x[[i]]$splits ))
    colnames(splits) <- paste0("split_", 1:ncol(splits))
  } else{
    splits <- NULL
  }

  # return instance of S3 class 'GenericML'
  return(
    structure(
      list(VEIN = vein,
           best = best.learners,
           propensity_scores = x[[1]]$propensity_scores,
           GenericML_single = NULL,
           splits = splits,
           generic_targets = generic_targets,
           arguments = list(learners_GenericML       = x[[1]]$arguments$learners_GenericML,
                            learner_propensity_score = x[[1]]$arguments$learner_propensity_score,
                            num_splits               = dim(generic_targets[[1]]$BLP)[3],
                            HT                       = x[[1]]$arguments$HT,
                            X1_BLP                   = x[[1]]$X1_BLP,
                            X1_GATES                 = x[[1]]$X1_GATES,
                            quantile_cutoffs         = x[[1]]$quantile_cutoffs,
                            diff_GATES               = x[[1]]$diff_GATES,
                            diff_CLAN                = x[[1]]$diff_CLAN,
                            vcov_BLP                 = x[[1]]$vcov_BLP,
                            vcov_GATES               = x[[1]]$vcov_GATES,
                            equal_variances_CLAN     = x[[1]]$equal_variances_CLAN,
                            prop_aux                 = x[[1]]$prop_aux,
                            significance_level       = x[[1]]$significance_level,
                            min_variation            = x[[1]]$min_variation,
                            parallel                 = NULL,
                            num_cores                = NULL,
                            seed                     = NULL,
                            store_learners           = NULL,
                            store_splits             = x[[1]]$store_splits)),
      class = "GenericML"))

} # FUN
