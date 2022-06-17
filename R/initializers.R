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


get_gates.3d <- function(num.learners, learners.names, num.generic_targets.gates){

  array(NA_real_, dim = c(num.generic_targets.gates, 7, num.learners),
        dimnames = list(paste0("gamma.", 1:num.generic_targets.gates), c("Estimate", "CB lower", "CB upper", "Std. Error", "z value", "Pr(<z)", "Pr(>z)"), learners.names))

} # FUN


get_best.3d <- function(num.learners, learners.names){

  array(NA_real_, dim = c(1, 2, num.learners),
        dimnames = list(NULL, c("lambda", "lambda.bar"), learners.names))

} # FUN


get_clan.3d.ls <- function(num.learners, learners.names, num.generic_targets.clan, num.vars.in.Z_CLAN, Z_CLAN.names){

  clan.3d.ls <- lapply(1:num.vars.in.Z_CLAN,
                       function(...) array(NA_real_, dim = c(num.generic_targets.clan, 7, num.learners), dimnames = list(paste0("gamma.", 1:num.generic_targets.clan), c("Estimate", "CB lower", "CB upper", "Std. Error", "z value", "Pr(<z)", "Pr(>z)"), learners.names)) )
  names(clan.3d.ls) <- Z_CLAN.names
  return(clan.3d.ls)

} # FUN


#' Setup function for \code{diff} arguments
#'
#' This setup function controls how differences of generic target parameters are taken. Returns a list with two components, called \code{subtract_from} and \code{subtracted}. The first element (\code{subtract_from}) denotes what shall be the base group to subtract from in the generic targets of interest (GATES or CLAN); either \code{"most"} or \code{"least"}. The second element (\code{subtracted}) are the groups to be subtracted from \code{subtract_from}, which is a subset of \eqn{{1,2,...,K}}, where \eqn{K} equals the number of groups. The number of groups should be consistent with the number of groups induced by the argument \code{quantile_cutoffs}, which is the cardinality of \code{quantile_cutoffs}, plus one.
#'
#' @details
#' The output of this setup function is intended to be used as argument in the functions \code{\link{GenericML}()} and \code{\link{GenericML_single}()} (arguments \code{diff_GATES}, \code{diff_CLAN}), as well as \code{\link{GATES}()} and \code{\link{CLAN}()} (argument \code{diff}).
#'
#' @param subtract_from String indicating the base group to subtract from, either \code{"most"} (default) or \code{"least"}. The most affected group corresponds to the \eqn{K}-th group in the paper (there are \eqn{K} groups). The least affected group corresponds to the first group.
#' @param subtracted Vector indicating the groups to be subtracted from the group specified in \code{subtract_from}. If there are \eqn{K} groups, \code{subtracted} should be a subset of \eqn{{1,2,...,K}}. Be careful to not specify a zero difference: If \code{subtract_from = "most"}, subtracting group K results in a zero difference. Same if \code{subtract_from = "least"} and we subtract group 1.
#'
#' @return
#' An object of class \code{"setup_diff"}, consisting of the following components:
#' \describe{
#'   \item{\code{subtract_from}}{A character equal to \code{"most"} or \code{"least"}.}
#'   \item{\code{subtracted}}{A numeric vector of group indices.}
#'   }
#' See the description above for details.
#'
#' @examples
#' ## specify quantile cutoffs (the 4 quartile groups here)
#' quantile_cutoffs <- c(0.25, 0.5, 0.75)
#'
#' ## Use group difference GK-G1 as generic targets in GATES and CLAN
#' ## Gx is the x-th group
#' setup_diff(subtract_from = "most", subtracted = 1)
#'
#' ## Use GK-G1, GK-G2, GK-G3 as differenced generic targets
#' setup_diff(subtract_from = "most", subtracted = c(1,2,3))
#'
#' ## Use G1-G2, G1-G3 as differenced generic targets
#' setup_diff(subtract_from = "least", subtracted = c(3,2))
#'
#' @references
#' Chernozhukov V., Demirer M., Duflo E., Fernández-Val I. (2020). \dQuote{Generic Machine Learning Inference on Heterogenous Treatment Effects in Randomized Experiments.} \emph{arXiv preprint arXiv:1712.04802}. URL: \url{https://arxiv.org/abs/1712.04802}.
#'
#' @seealso
#' \code{\link{GenericML}()},
#' \code{\link{GenericML_single}()},
#' \code{\link{CLAN}()},
#' \code{\link{GATES}()},
#' \code{\link{setup_X1}()},
#' \code{\link{setup_vcov}()}
#'
#' @export
setup_diff <- function(subtract_from = "most",
                       subtracted = 1){

  stopifnot(subtract_from %in% c("most", "least") & length(subtract_from) == 1)
  stopifnot(is.vector(subtracted) & is.numeric(subtracted))

  structure(
    list(subtract_from = subtract_from,
         subtracted = subtracted),
    class = "setup_diff")

} # FUN


#' Setup function for \code{vcov_control} arguments
#'
#' Returns a list with two elements called \code{estimator} and \code{arguments}. The element \code{estimator} is a string specifying the covariance matrix estimator to be used in the linear regression regression of interest and needs to be a covariance estimator function in the \href{https://CRAN.R-project.org/package=sandwich}{"sandwich"} package. The second element, \code{arguments}, is a list of arguments that shall be passed to the function specified in the first element, \code{estimator}.
#'
#' @details
#' The output of this setup function is intended to be used as argument in the functions \code{\link{GenericML}()} and \code{\link{GenericML_single}()} (arguments \code{vcov_BLP}, \code{vcov_GATES}), as well as \code{\link{BLP}()} and \code{\link{GATES}()} (argument \code{vcov_control}).
#'
#' @param estimator Character specifying a covariance matrix estimator in the \href{https://CRAN.R-project.org/package=sandwich}{"sandwich"} package. Default is \code{"vcovHC"}. Supported estimators are \code{"vcovBS"}, \code{"vcovCL"}, \code{"vcovHAC"}, and \code{"vcovHC"}.
#' @param arguments A list of arguments that are to be passed to the function in the \code{"sandwich"} package that is specified in \code{estimator}. Default is \code{list(type = "const")}, which specifies the homoskedastic ordinary least squares covariance matrix estimator.
#'
#' @return
#' An object of class \code{"setup_vcov"}, consisting of the following components:
#' \describe{
#'   \item{\code{estimator}}{A character equal to covariance estimation function names in the \href{https://CRAN.R-project.org/package=sandwich}{"sandwich"} package.}
#'   \item{\code{arguments}}{A list of arguments that shall be passed to the function specified in the \code{estimator} argument.}
#'   }
#' See the description above for details.
#'
#' @references
#' Zeileis A. (2004). \dQuote{Econometric Computing with HC and HAC Covariance Matrix Estimators.} \emph{Journal of Statistical Software}, \bold{11}(10), 1--17. \doi{10.18637/jss.v011.i10}
#'
#' Zeileis A. (2006). \dQuote{Object-Oriented Computation of Sandwich Estimators.} \emph{Journal of Statistical Software}, \bold{16}(9), 1--16. \doi{10.18637/jss.v016.i09}
#'
#' @seealso
#' \code{\link{GenericML}()},
#' \code{\link{GenericML_single}()},
#' \code{\link{BLP}()},
#' \code{\link{GATES}()},
#' \code{\link{setup_X1}()},
#' \code{\link{setup_diff}()}
#'
#' @examples
#' # use standard homoskedastic OLS covariance matrix estimate
#' setup_vcov(estimator = "vcovHC", arguments = list(type = "const"))
#'
#' # use White's heteroskedasticity-robust estimator
#' setup_vcov(estimator = "vcovHC", arguments = list(type = "HC0"))
#'
#' if (require("sandwich")){
#'
#' # use HAC-robust estimator with prewhitening and Andrews' (Econometrica, 1991) weights
#' # since weightsAndrews() is a function in 'sandwich', require this package
#' setup_vcov(estimator = "vcovHAC", arguments = list(prewhite = TRUE, weights = weightsAndrews))
#'
#' }
#'
#' @import sandwich
#'
#' @export
setup_vcov <- function(estimator = "vcovHC",
                       arguments = list(type = "const")){

  # input checks
  stopifnot(is.character(estimator))
  stopifnot(length(estimator) == 1)
  stopifnot(is.list(arguments))

  # check if optional arguments that will require subsetting later are of correct data type
  args_nam <- names(arguments)

  if("cluster" %in% args_nam){
    stopifnot(is.matrix(arguments$cluster) | is.vector(arguments$cluster))
    stopifnot(is.numeric(arguments$cluster))
  } else if("order.by" %in% args_nam){
    stopifnot(is.numeric(arguments$order.by) & is.vector(arguments$order.by))
  } else if("omega" %in% args_nam){
    stopifnot(is.numeric(arguments$omega) & is.vector(arguments$omega))
  } # IF


  if(!estimator %in% c("vcovBS", "vcovCL", "vcovHAC", "vcovHC")){
    stop(paste0("'estimator' needs to be in ",
         "c('vcovBS', 'vcovCL', 'vcovHAC', 'vcovHC')"), call. = FALSE)
  } # IF


  # return
  setup_vcov_NoChecks(
    estimator = estimator,
    arguments = arguments)

} # FUN


#' same as above, just without input checks
#'
#' @import sandwich
#' @noRd
setup_vcov_NoChecks <- function(estimator = "vcovHC",
                                arguments = list(type = "const")){
  structure(
    list(estimator = estimator,
         arguments = arguments),
    class = "setup_vcov")

} # FUN


# helper function that makes a 'setup_vcov' ready for further use
# x is of class "setup_vcov"
setup_vcov_align <- function(x){

  # this function simply makes the cluster argument a matrix (important for later use)
  if(!is.null(x$arguments$cluster)){
    x$arguments$cluster <- as.matrix(x$arguments$cluster)
  } # IF

  # recover class
  structure(x, class = "setup_vcov")

} # FUN


# helper function that subsets optional arguments that were passed in setup_vcov()
# x is of class "setup_vcov"
setup_vcov_subset <- function(x, idx){

  # subset the (potential) inputs that need subsetting due to sample splits
  # if an input is NULL, then the corresponding line won't have an effect (naturally)
  x$arguments$cluster  <- x$arguments$cluster[idx,,drop = FALSE] # matrix, rest are vectors
  x$arguments$order.by <- x$arguments$order.by[idx]
  x$arguments$omega    <- x$arguments$omega[idx]

  # recover class
  structure(x, class = "setup_vcov")

} # FUN


#' Setup function controlling the matrix \eqn{X_1} in the BLP or GATES regression
#'
#' Returns a list with three elements. The first element of the list, \code{funs_Z}, controls which functions of matrix \code{Z} are used as regressors in \eqn{X_1}. The second element, \code{covariates}, is an optional matrix of custom covariates that shall be included in \eqn{X_1}. The third element, \code{fixed_effects}, controls the inclusion of fixed effects.
#'
#' @details
#' The output of this setup function is intended to be used as argument in the functions \code{\link{GenericML}()} and \code{\link{GenericML_single}()} (arguments \code{X1_BLP}, \code{X1_GATES}), as well as \code{\link{BLP}()} and \code{\link{GATES}()} (argument \code{X1_control}).
#'
#' @param funs_Z Character vector controlling the functions of \code{Z} to be included in \eqn{X_1}. Subset of \code{c("S", "B", "p")}, where \code{"p"} corresponds to the propensity scores, \code{"B"} to the proxy baseline estimates, and \code{"S"} to the proxy CATE estimates. Default is \code{"B"}.
#' @param covariates Optional numeric matrix containing additional covariates to be included in \eqn{X_1}. Default is \code{NULL}.
#' @param fixed_effects Numeric vector of integers that indicates cluster membership of the observations: For each cluster, a fixed effect will be added. Default is \code{NULL} for no fixed effects.
#'
#' @return
#' An object of class \code{"setup_X1"}, consisting of the following components:
#' \describe{
#'   \item{\code{funs_Z}}{A character vector, being a subset of \code{c("S", "B", "p")}.}
#'   \item{\code{covariates}}{Either \code{NULL} or a numeric matrix.}
#'   \item{\code{fixed_effects}}{Either \code{NULL} or an integer vector indicating cluster membership.}
#'   }
#' See the description above for details.
#'
#' @references
#' Chernozhukov V., Demirer M., Duflo E., Fernández-Val I. (2020). \dQuote{Generic Machine Learning Inference on Heterogenous Treatment Effects in Randomized Experiments.} \emph{arXiv preprint arXiv:1712.04802}. URL: \url{https://arxiv.org/abs/1712.04802}.
#'
#' @seealso
#' \code{\link{GenericML}()},
#' \code{\link{GenericML_single}()},
#' \code{\link{BLP}()},
#' \code{\link{GATES}()},
#' \code{\link{setup_vcov}()},
#' \code{\link{setup_diff}()}
#'
#' @examples
#' set.seed(1)
#' n <- 100 # sample size
#' p <- 5   # number of covariates
#' covariates <- matrix(runif(n*p), n, p) # sample matrix of covariates
#'
#' # let there be three clusters; assign membership randomly
#' fixed_effects <- sample(c(1,2,3), size = n, replace = TRUE)
#'
#' # use BCA estimates in matrix X1
#' setup_X1(funs_Z = "B", covariates = NULL, fixed_effects = NULL)
#'
#' # use BCA and propensity score estimates in matrix X1
#' # uses uniform covariates and fixed effects
#' setup_X1(funs_Z = c("B", "p"), covariates = covariates, fixed_effects = NULL)
#'
#' @export
setup_X1 <- function(funs_Z = c("B"),
                     covariates = NULL,
                     fixed_effects = NULL){

  # input checks
  legalinput <- funs_Z %in% c("S", "B", "p")

  if(!all(legalinput)){

    stop(paste0("Entries '",
                paste(funs_Z[!legalinput], collapse = "', '"),
                "' of ",
                "funs_Z are not contained in c('S', 'B', 'p').",
                " The entries must be a subset of this vector."), call. = FALSE)

  } # IF


  if(!is.null(covariates)){

    if(!(is.numeric(covariates) & is.matrix(covariates))){
      stop("If supplied, 'covariates' must be a numeric matrix. Did you supply a data frame?",
           call. = FALSE)
    } # IF

    if(any(is.na(covariates))) stop("'covariates' contains missing values", call. = FALSE)

  } # IF !NULL


  if(!is.null(fixed_effects)){

    if(!(is.numeric(fixed_effects) & is.vector(fixed_effects))){
      stop("If supplied, 'fixed_effects' must be a numeric vector",
           call. = FALSE)
    } # IF

    if(any(fixed_effects %% 1 != 0)){
      stop("All elements in the vector 'fixed_effects' must be integer-valued",
           call. = FALSE)
    } # IF


  } # IF !NULL


  # return
  setup_X1_NoChecks(funs_Z = funs_Z,
                    covariates = covariates,
                    fixed_effects = fixed_effects)

} # FUN


# same as above, just without input checks
setup_X1_NoChecks <- function(funs_Z = c("B"),
                              covariates = NULL,
                              fixed_effects = NULL){

  if(is.null(fixed_effects)){
    temp <- NULL
  } else{
    temp <- as.integer(fixed_effects)
  }

  structure(
    list(funs_Z = funs_Z,
         covariates = covariates,
         fixed_effects = temp),
    class = "setup_X1")

} # FUN
