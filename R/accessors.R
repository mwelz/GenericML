#' Accessor function for the BLP generic target estimates
#'
#' @param x An object of the class \code{"\link{GenericML}"}, as returned by the function \code{\link{GenericML}()}.
#' @param learner A character string of the learner whose BLP generic target estimates shall be accessed. Default is \code{"best"} for the best learner for BLP.
#' @param plot Logical. If \code{TRUE} (default), a \code{"\link[ggplot2]{ggplot}"} object is computed.
#'
#' @return An object of class \code{"BLP_info"}, which consists of the following components:
#' \describe{
#'   \item{\code{estimate}}{A numeric vector of point estimates of the BLP generic targets.}
#'   \item{\code{confidence_interval}}{A numeric matrix of the lower and upper confidence bounds for each generic target. The confidence level of the implied confidence interval is equal to \code{1 - 2 * significance_level}.}
#'   \item{\code{confidence_level}}{The confidence level of the confidence intervals. Equals \code{1 - 2 * significance_level}.}
#'   \item{\code{learner}}{The argument \code{learner}.}
#'   \item{\code{plot}}{An object of class \code{"\link[ggplot2]{ggplot}"}. Only returned if the argument \code{plot = TRUE}.}
#'   }
#'
#' @examples
#' if(require("rpart") && require("ranger")){
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
#' learners <- c("tree", "mlr3::lrn('ranger', num.trees = 10)")
#'
#' ## perform generic ML inference
#' # small number of splits to keep computation time low
#' x <- GenericML(Z, D, Y, learners, num_splits = 2,
#'                parallel = FALSE)
#'
#' ## access best learner
#' get_best(x)
#'
#' ## access BLP generic targets for best learner w/o plot
#' get_BLP(x, learner = "best", plot = FALSE)
#'
#' ## access BLP generic targets for ranger learner w/o plot
#' get_BLP(x, learner = "mlr3::lrn('ranger', num.trees = 10)", plot = FALSE)
#'
#' ## access GATES generic targets for best learner w/o plot
#' get_GATES(x, learner = "best", plot = FALSE)
#'
#' ## access GATES generic targets for ranger learner w/o plot
#' get_GATES(x, learner = "mlr3::lrn('ranger', num.trees = 10)", plot = FALSE)
#'
#' ## access CLAN generic targets for "V1" & best learner, w/o plot
#' get_CLAN(x, learner = "best", variable = "V1", plot = FALSE)
#'
#' ## access CLAN generic targets for "V1" & ranger learner, w/o plot
#' get_CLAN(x, learner = "mlr3::lrn('ranger', num.trees = 10)",
#'          variable = "V1", plot = FALSE)
#' }
#'
#' @seealso
#' \code{\link{GenericML}()},
#' \code{\link{get_GATES}()},
#' \code{\link{get_CLAN}()},
#' \code{\link{get_best}()},
#' \code{\link{print.BLP_info}()},
#' \code{\link{print.GATES_info}()},
#' \code{\link{print.CLAN_info}()}
#'
#' @export
get_BLP <- function(x, learner = "best", plot = TRUE){

  # access the results
  results <- accessor_BLP_GATES(x = x, type = "BLP", learner = learner)

  # prepare the output
  out <- accessor_output(x = x, accessor_obj = results, plot = plot,
                         type = "BLP", learner = learner,
                         CLAN_variable = NULL, ATE = TRUE)

  # assign class and return
  class(out) <- "BLP_info"
  return(out)

} # FUN




#' Accessor function for the GATES generic target estimates
#'
#' @param x An object of the class \code{"\link{GenericML}"}, as returned by the function \code{\link{GenericML}()}.
#' @param learner A character string of the learner whose GATES generic target estimates shall be accessed. Default is \code{"best"} for the best learner for GATES.
#' @param plot Logical. If \code{TRUE} (default), a \code{"\link[ggplot2]{ggplot}"} object is computed.
#'
#' @return An object of class \code{"GATES_info"}, which consists of the following components:
#' \describe{
#'   \item{\code{estimate}}{A numeric vector of point estimates of the GATES generic targets.}
#'   \item{\code{confidence_interval}}{A numeric matrix of the lower and upper confidence bounds for each generic target. The confidence level of the implied confidence interval is equal to \code{1 - 2 * significance_level}.}
#'   \item{\code{confidence_level}}{The confidence level of the confidence intervals. Equals \code{1 - 2 * significance_level}.}
#'   \item{\code{learner}}{The argument \code{learner}.}
#'   \item{\code{plot}}{An object of class \code{"\link[ggplot2]{ggplot}"}. Only returned if the argument \code{plot = TRUE}.}
#'   }
#'
#' @examples
#' if(require("rpart") && require("ranger")){
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
#' learners <- c("tree", "mlr3::lrn('ranger', num.trees = 10)")
#'
#' ## perform generic ML inference
#' # small number of splits to keep computation time low
#' x <- GenericML(Z, D, Y, learners, num_splits = 2,
#'                parallel = FALSE)
#'
#' ## access best learner
#' get_best(x)
#'
#' ## access BLP generic targets for best learner w/o plot
#' get_BLP(x, learner = "best", plot = FALSE)
#'
#' ## access BLP generic targets for ranger learner w/o plot
#' get_BLP(x, learner = "mlr3::lrn('ranger', num.trees = 10)", plot = FALSE)
#'
#' ## access GATES generic targets for best learner w/o plot
#' get_GATES(x, learner = "best", plot = FALSE)
#'
#' ## access GATES generic targets for ranger learner w/o plot
#' get_GATES(x, learner = "mlr3::lrn('ranger', num.trees = 10)", plot = FALSE)
#'
#' ## access CLAN generic targets for "V1" & best learner, w/o plot
#' get_CLAN(x, learner = "best", variable = "V1", plot = FALSE)
#'
#' ## access CLAN generic targets for "V1" & ranger learner, w/o plot
#' get_CLAN(x, learner = "mlr3::lrn('ranger', num.trees = 10)",
#'          variable = "V1", plot = FALSE)
#' }
#'
#' @seealso
#' \code{\link{GenericML}()},
#' \code{\link{get_BLP}()},
#' \code{\link{get_CLAN}()},
#' \code{\link{get_best}()},
#' \code{\link{print.BLP_info}()},
#' \code{\link{print.GATES_info}()},
#' \code{\link{print.CLAN_info}()}
#'
#' @export
get_GATES <- function(x, learner = "best", plot = TRUE){

  # access the results
  results <- accessor_BLP_GATES(x = x, type = "GATES", learner = learner)

  # prepare the output
  out <- accessor_output(x = x, accessor_obj = results, plot = plot,
                         type = "GATES", learner = learner,
                         CLAN_variable = NULL, ATE = TRUE)

  # assign class and return
  class(out) <- "GATES_info"
  return(out)

} # FUN


#' Accessor function for the CLAN generic target estimates
#'
#' @param x An object of the class \code{"\link{GenericML}"}, as returned by the function \code{\link{GenericML}()}.
#' @param variable The (character) name of a variabe on which CLAN was performed.
#' @param learner A character string of the learner whose CLAN generic target estimates shall be accessed. Default is \code{"best"} for the best learner for CLAN
#' @param plot Logical. If \code{TRUE} (default), a \code{"\link[ggplot2]{ggplot}"} object is computed.
#'
#' @return An object of class \code{"CLAN_info"}, which consists of the following components:
#' \describe{
#'   \item{\code{estimate}}{A numeric vector of point estimates of the CLAN generic targets.}
#'   \item{\code{confidence_interval}}{A numeric matrix of the lower and upper confidence bounds for each generic target. The confidence level of the implied confidence interval is equal to \code{1 - 2 * significance_level}.}
#'   \item{\code{confidence_level}}{The confidence level of the confidence intervals. Equals \code{1 - 2 * significance_level}.}
#'   \item{\code{learner}}{The argument \code{learner}.}
#'   \item{\code{plot}}{An object of class \code{"\link[ggplot2]{ggplot}"}. Only returned if the argument \code{plot = TRUE}.}
#'   \item{\code{CLAN_variable}}{The name of the CLAN variable of interest.}
#'   }
#'
#' @examples
#' if(require("rpart") && require("ranger")){
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
#' learners <- c("tree", "mlr3::lrn('ranger', num.trees = 10)")
#'
#' ## perform generic ML inference
#' # small number of splits to keep computation time low
#' x <- GenericML(Z, D, Y, learners, num_splits = 2,
#'                parallel = FALSE)
#'
#' ## access best learner
#' get_best(x)
#'
#' ## access BLP generic targets for best learner w/o plot
#' get_BLP(x, learner = "best", plot = FALSE)
#'
#' ## access BLP generic targets for ranger learner w/o plot
#' get_BLP(x, learner = "mlr3::lrn('ranger', num.trees = 10)", plot = FALSE)
#'
#' ## access GATES generic targets for best learner w/o plot
#' get_GATES(x, learner = "best", plot = FALSE)
#'
#' ## access GATES generic targets for ranger learner w/o plot
#' get_GATES(x, learner = "mlr3::lrn('ranger', num.trees = 10)", plot = FALSE)
#'
#' ## access CLAN generic targets for "V1" & best learner, w/o plot
#' get_CLAN(x, learner = "best", variable = "V1", plot = FALSE)
#'
#' ## access CLAN generic targets for "V1" & ranger learner, w/o plot
#' get_CLAN(x, learner = "mlr3::lrn('ranger', num.trees = 10)",
#'          variable = "V1", plot = FALSE)
#' }
#'
#' @seealso
#' \code{\link{GenericML}()},
#' \code{\link{get_BLP}()},
#' \code{\link{get_GATES}()},
#' \code{\link{get_best}()},
#' \code{\link{print.BLP_info}()},
#' \code{\link{print.GATES_info}()},
#' \code{\link{print.CLAN_info}()}
#'
#' @export
get_CLAN <- function(x, variable, learner = "best", plot = TRUE){

  # access the results
  results <- accessor_CLAN(x = x, variable = variable, learner = learner)

  # prepare the output
  out <- accessor_output(x = x, accessor_obj = results, plot = plot,
                         type = "CLAN", learner = learner,
                         CLAN_variable = variable, ATE = FALSE)

  # assign class and return
  class(out) <- "CLAN_info"
  return(out)

} # FUN


#' An internal accessor function for CLAN analyses
#'
#' @param x A \code{"GenericML"} object
#' @param variable Name of variable for which CLAN should be performed
#' @param learner Learner of the analysis, either \code{"best"} or learners used in \code{x} (error if not)
#'
#' @return A matrix of point estimates, confidence intervals, and adjusted p values (upper and lower)
#' @noRd
accessor_CLAN <- function(x, variable, learner)
{
  # input check
  isGenericMLcheck(x)

  # get CLAN names
  CLAN_vars <- CLAN_names(x)

  if(!(variable %in% CLAN_vars)){
    stop(paste0("No CLAN was performed on this variable. ",
                "CLAN was performed on the variables ",
                paste0(CLAN_vars, collapse = ", ")))
  } # IF

  # call the main function
  accessor_CLAN_noChecks(x = x,
                         variable = variable,
                         learner = learner)

} # FUN


# same as accessor_CLAN(), just without input checks
accessor_CLAN_noChecks <- function(x, variable, learner)
{
  if(identical(learner, "best")){

    out <- x$VEIN$best_learners$CLAN[[variable]]

  } else{

    if(!(learner %in% x$arguments$learners_GenericML)){

      stop("Specified learner is not used in this instance of GenericML")

    } # IF

    out <- x$VEIN$all_learners$CLAN[[learner]][[variable]]

  } # IF

  return(out)
} # FUN


#' An internal accessor function for BLP or GATES analyses
#'
#' @param x A \code{"GenericML"} object
#' @param type Type of analysis, either \code{"BLP"} or \code{"GATES"}
#' @param learner Learner of the analysis, either \code{"best"} or learners used in \code{x} (error if not)
#'
#' @return A matrix of point estimates, confidence intervals, and adjusted p values (upper and lower)
#' @noRd
accessor_BLP_GATES <- function(x, type, learner)
{
  # input check
  isGenericMLcheck(x)

  if(identical(learner, "best")){

    out <- x$VEIN$best_learners[[type]]

  } else{

    if(!(learner %in% x$arguments$learners_GenericML)){

      stop("Specified learner is not used in this instance of GenericML")

    } # IF

    out <- x$VEIN$all_learners[[type]][[learner]]

  } # IF

  return(out)

} # FUN


#' Internal function that prepares the output of accessor function
#'
#' @param x An object of class \code{"GenericML"}
#' @param accessor_obj An object as returned by \code{accessor_CLAN()} or \code{accessor_BLP_GATES()}
#' @param plot Logical; shall \code{"ggplot"} object be included in the output?
#' @param type Either \code{"BLP"}, \code{"GATES"}, or \code{"CLAN"}
#' @param learner Learner of the analysis, possibly \code{"best"}
#' @param CLAN_variable Variable along which CLAN shall be performed. Only applicable if \code{type = "CLAN"}
#' @param ATE Shall ATE be included in plot?
#'
#' @return A list of point estimates, confidence intervals, and adjusted p values (minimum of the lower and upper p value estimates)
#' @noRd
accessor_output <- function(x, accessor_obj, plot, type, learner, CLAN_variable, ATE)
{
  ## NOTE: the checks have already been performed in accessor_CLAN() or accessor_BLP_GATES()
  ## prepare components to return
  Estimate <- accessor_obj[, "Estimate"]
  ConfidenceInterval <- accessor_obj[, c("CB lower", "CB upper")]
  colnames(ConfidenceInterval) <- c("lower", "upper")

  ## only return the minimum of the two probabilities as p value
  pval <- pmin(accessor_obj[, "Pr(<z) adjusted"], accessor_obj[, "Pr(>z) adjusted"])

  ## prepare output
  # Generic ML estimation has size control of 2 * significance_level
  out <- list(estimate = Estimate, confidence_interval = ConfidenceInterval, p_value = pval,
              confidence_level = confidence_level(x), learner = learner)

  ## if requested, add ggplot object to output
  if(isTRUE(plot)){
    out$plot <- plot.GenericML(x = x, type = type,
                               learner = learner, CLAN_variable = CLAN_variable, ATE = ATE)
  } # IF

  ## if requested, add CLAN variable to output
  if(!is.null(CLAN_variable)){
    out$CLAN_variable <- CLAN_variable
  } # IF

  return(out)

} # FUN


#' Accessor function for the best learner estimates
#'
#' The best learner is determined by maximizing the criteria \eqn{\Lambda} and \eqn{\bar{\Lambda}}, see Sections 5.2 and 5.3 of the paper. This function accesses the estimates of these two criteria,
#'
#' @param x An object of the class \code{"\link{GenericML}"}, as returned by the function \code{\link{GenericML}()}.
#'
#' @return
#' An object of class \code{"best"}, which consists of the following components:
#' \describe{
#'   \item{\code{BLP}}{A string holding the name of the best learner for a BLP analysis.}
#'   \item{\code{GATES}}{A string holding the name of the best learner for a GATES analysis.}
#'   \item{\code{CLAN}}{A string holding the name of the best learner for a CLAN analysis (same learner as in \code{GATES}).}
#'   \item{\code{overview}}{A numeric matrix of the estimates of the performance measures \eqn{\Lambda} and \eqn{\bar{\Lambda}} for each learner.}}
#'
#' @examples
#' if(require("rpart") && require("ranger")){
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
#' learners <- c("tree", "mlr3::lrn('ranger', num.trees = 10)")
#'
#' ## perform generic ML inference
#' # small number of splits to keep computation time low
#' x <- GenericML(Z, D, Y, learners, num_splits = 2,
#'                parallel = FALSE)
#'
#' ## access best learner
#' get_best(x)
#'
#' ## access BLP generic targets for best learner w/o plot
#' get_BLP(x, learner = "best", plot = FALSE)
#'
#' ## access BLP generic targets for ranger learner w/o plot
#' get_BLP(x, learner = "mlr3::lrn('ranger', num.trees = 10)", plot = FALSE)
#'
#' ## access GATES generic targets for best learner w/o plot
#' get_GATES(x, learner = "best", plot = FALSE)
#'
#' ## access GATES generic targets for ranger learner w/o plot
#' get_GATES(x, learner = "mlr3::lrn('ranger', num.trees = 10)", plot = FALSE)
#'
#' ## access CLAN generic targets for "V1" & best learner, w/o plot
#' get_CLAN(x, learner = "best", variable = "V1", plot = FALSE)
#'
#' ## access CLAN generic targets for "V1" & ranger learner, w/o plot
#' get_CLAN(x, learner = "mlr3::lrn('ranger', num.trees = 10)",
#'          variable = "V1", plot = FALSE)
#' }
#'
#' @seealso
#' \code{\link{GenericML}()},
#' \code{\link{get_BLP}()},
#' \code{\link{get_GATES}()},
#' \code{\link{get_CLAN}()}
#'
#' @export
get_best <- function(x)
{
  isGenericMLcheck(x)
  out <- x$best
  out$num_splits <- x$arguments$num_splits
  return(structure(out, class = "best"))
} # FUN
