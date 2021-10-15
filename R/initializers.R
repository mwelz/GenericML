initializer.for.splits <- function(Z, Z_CLAN, learners,
                                   num_splits, quantile_cutoffs,
                                   diff_GATES,
                                   diff_CLAN){

  # helper function that initializes object in the generic ML splitting procedure

  K <- length(quantile_cutoffs) + 1
  CLAN_group.base  <- ifelse(diff_CLAN$subtract_from == "least", 1, K)
  GATES_group.base <- ifelse(diff_GATES$subtract_from == "least", 1, K)
  CLAN_subtracted  <- diff_CLAN$subtracted
  GATES_subtracted <- diff_GATES$subtracted

  clan <- array(NA_real_, dim = c(K + length(CLAN_subtracted), 7, num_splits),
                dimnames = list(c(paste0("delta.", 1:K),
                                  paste0(
                                    "delta.", CLAN_group.base, "-",
                                    "delta.", CLAN_subtracted)),
                                c("Estimate", "CB lower", "CB upper", "Std. Error",
                                  "z value", "Pr(<z)", "Pr(>z)"),
                                NULL))

  gates <- array(NA_real_, dim = c(K + length(GATES_subtracted), 7, num_splits),
                 dimnames = list(c(paste0("gamma.", 1:K),
                                   paste0(
                                     "gamma.", GATES_group.base, "-",
                                     "gamma.", GATES_subtracted)),
                                 c("Estimate", "CB lower", "CB upper", "Std. Error",
                                   "z value", "Pr(<z)", "Pr(>z)"),
                                 NULL))

  blp <- array(NA_real_, dim = c(2, 7, num_splits),
               dimnames = list(c("beta.1", "beta.2"),
                               c("Estimate", "CB lower", "CB upper", "Std. Error",
                                 "z value", "Pr(<z)", "Pr(>z)"),
                               NULL))

  best <- array(NA_real_, dim = c(1, 2, num_splits),
                dimnames = list(NULL, c("lambda", "lambda.bar"), NULL))

  clan.lists <- lapply(1:ncol(Z_CLAN), function(...) clan )
  names(clan.lists) <- colnames(Z_CLAN)

  out.ls <- lapply(1:length(learners), function(...){

    list(BLP = blp, GATES = gates, CLAN = clan.lists, best = best)

  })

  names(out.ls) <- learners
  return(out.ls)

} # END FUN



initialize.gen.ml <- function(generic.ml.across.learners.obj){

  # helper function for initialization of the final returned object

  gates.nam  <- rownames(generic.ml.across.learners.obj[[1]]$GATES[,,1])
  z.clan.nam <- names(generic.ml.across.learners.obj[[1]]$CLAN)
  clan.nam   <- rownames(generic.ml.across.learners.obj[[1]]$CLAN[[1]][,,1])
  learners   <- names(generic.ml.across.learners.obj)

  blp.mat <- matrix(NA_real_, nrow = 2, ncol = 5,
                    dimnames = list(c("beta.1", "beta.2"),
                                    c("Estimate", "CB lower", "CB upper",
                                      "Pr(<z) adjusted", "Pr(>z) adjusted")))

  gates.mat <- matrix(NA_real_, nrow = length(gates.nam), ncol = 5,
                      dimnames = list(gates.nam,
                                      c("Estimate", "CB lower", "CB upper",
                                        "Pr(<z) adjusted", "Pr(>z) adjusted")))

  clan.mat <- matrix(NA_real_, nrow = length(clan.nam), ncol = 5,
                     dimnames = list(clan.nam,
                                     c("Estimate", "CB lower", "CB upper",
                                       "Pr(<z) adjusted", "Pr(>z) adjusted")))


  gates.ls <- lapply(learners, function(...) gates.mat)
  names(gates.ls) <- learners

  blp.ls <- lapply(learners, function(...) blp.mat)
  names(blp.ls) <- learners

  z.clan.ls <- lapply(z.clan.nam, function(...) clan.mat)
  names(z.clan.ls) <- z.clan.nam

  clan.ls <- lapply(learners, function(...) z.clan.ls)
  names(clan.ls) <- learners

  return(list(BLP = blp.ls, GATES = gates.ls, CLAN = clan.ls))

} # END FUN


# for parallelized across learner function:
get_blp.3d <- function(num.learners, learners.names){

  array(NA_real_, dim = c(2, 7, num.learners), dimnames = list(c("beta.1", "beta.2"), c("Estimate", "CB lower", "CB upper", "Std. Error", "z value", "Pr(<z)", "Pr(>z)"), learners.names))

} # FUN


get_gates.3d <- function(num.learners, learners.names, num.generic.targets.gates){

  array(NA_real_, dim = c(num.generic.targets.gates, 7, num.learners),
        dimnames = list(paste0("gamma.", 1:num.generic.targets.gates), c("Estimate", "CB lower", "CB upper", "Std. Error", "z value", "Pr(<z)", "Pr(>z)"), learners.names))

} # FUN


get_best.3d <- function(num.learners, learners.names){

  array(NA_real_, dim = c(1, 2, num.learners),
        dimnames = list(NULL, c("lambda", "lambda.bar"), learners.names))

} # FUN


get_clan.3d.ls <- function(num.learners, learners.names, num.generic.targets.clan, num.vars.in.Z_CLAN, Z_CLAN.names){

  clan.3d.ls <- lapply(1:num.vars.in.Z_CLAN,
                       function(...) array(NA_real_, dim = c(num.generic.targets.clan, 7, num.learners), dimnames = list(paste0("gamma.", 1:num.generic.targets.clan), c("Estimate", "CB lower", "CB upper", "Std. Error", "z value", "Pr(<z)", "Pr(>z)"), learners.names)) )
  names(clan.3d.ls) <- Z_CLAN.names
  return(clan.3d.ls)

} # FUN


#' Setup function for the arguments \code{diff_GATES} and \code{diff_CLAN} in the \code{\link{GenericML}} function.
#'
#' Returns a list with two elements called \code{subtract_from} and \code{subtracted}. The first element (\code{subtract_from}) denotes what shall be the base group to subtract from in the generic targets of interest (GATES or CLAN); either \code{"most"} or \code{"least"}. The second element (\code{subtracted}) are the groups to be subtracted from \code{subtract_from}, which is a subset of \eqn{{1,2,...,K}}, where \eqn{K} equals the number of groups. The number of groups should be consistent with the number of groups induced by the argument \code{quantile_cutoffs}.
#'
#' @param subtract_from String indicating the base group to subtract from, either \code{"most"} (default) or \code{"least"}. The most affected group corresponds to the \eqn{K}-th group in the paper (there are \eqn{K} groups). The least affected group corresponds to the first group.
#' @param subtracted Vector indicating the groups to be subtracted from the group specified in \code{subtract_from}. If there are \eqn{K} group, \code{subtracted} should be a subset of \eqn{{1,2,...,K}}.
#'
#' @export
setup_diff <- function(subtract_from = "most",
                            subtracted = 1){

  list(subtract_from = subtract_from,
       subtracted = subtracted)

} # FUN


#' Setup function for the arguments \code{vcov_BLP} and \code{vcov_GATES} in the \code{\link{GenericML}} function.
#'
#' Returns a list with two elements called \code{estimator} and \code{arguments}. The element \code{estimator} is a string specifying the covariance matrix estimator to be used in the linear regression regression of interest and needs to be a covariance estimator function in the \href{https://cran.r-project.org/web/packages/sandwich/}{sandwich} package. The second element, \code{arguments}, is a list of arguments that shall be passed to the function specified in the first element, \code{estimator}.
#'
#' @param estimator String specifying a covariance matrix estimator in the \href{https://cran.r-project.org/web/packages/sandwich/}{sandwich}. Default is \code{"vcovHC"}. Recommended estimators are \code{"vcovBS"}, \code{"vcovCL"}, \code{"vcovHAC"}, and \code{"vcovHC"}.
#' @param arguments A list of arguments that are to be passed to the function in the sandwich package that is specified in \code{estimator}. Default is \code{list(type = "const")}, which specifies the homoskedastic ordinary least squares covariance matrix estimate.
#'
#' @export
setup_vcov <- function(estimator = "vcovHC",
                       arguments = list(type = "const")){

  list(estimator = estimator,
       arguments = arguments)

} # FUN
