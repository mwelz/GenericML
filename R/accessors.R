#' Accessor function for the BLP generic target estimates
#'
#' @param x An object of the class \code{\link{GenericML}}.
#' @param learner A character string of the learner whose BLP generic target estimates shall be accessed. Default is \code{"best"} for the best learner for BLP.
#' @param plot Logical. If \code{TRUE} (default), a plot is printed.
#'
#' @return A numeric matrix of BLP generic target estimates which contains information on point estimates, confidence bounds, and (adjusted) p-values. Furthermore, prints a plot if \code{plot = TRUE}.
#'
#' @examples
#' if(require("glmnet") && require("ranger")){
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
#' \code{\link{GenericML}},
#' \code{\link{get_GATES}},
#' \code{\link{get_CLAN}},
#' \code{\link[=plot.GenericML]{plot}}
#'
#' @export
get_BLP <- function(x, learner = "best", plot = TRUE){

  if(class(x) != "GenericML") stop("x needs to be instance of the class GenericML")

  if(learner == "best"){

    temp <- x$VEIN$best_learners$BLP

  } else{

    if(!(learner %in% x$arguments$learners_GenericML)){

      stop("Specified learner is not used in this instance of GenericML")

    } # IF

    temp <- x$VEIN$all_learners$BLP[[learner]]

  } # IF

  # only print the minimum of the two probabilities as p-value
  pval <- pmin(temp[, "Pr(<z) adjusted"], temp[, "Pr(>z) adjusted"])
  out  <- cbind(temp[, c("Estimate", "CB lower", "CB upper")],
                pval)
  colnames(out)[4] <- "Pr(>|z|)"


  if(plot) print(plot.GenericML(x = x, learner = learner, type = "BLP"))

  out

} # FUN




#' Accessor function for the GATES generic target estimates
#'
#' @param x An object of the class \code{\link{GenericML}}.
#' @param learner A character string of the learner whose GATES generic target estimates shall be accessed. Default is \code{"best"} for the best learner for GATES.
#' @param plot Logical. If \code{TRUE} (default), a plot is printed.
#'
#' @return A numeric matrix of GATES generic target estimates which contains information on point estimates, confidence bounds, and (adjusted) p-values. Furthermore, prints a plot if \code{plot = TRUE}.
#'
#' @examples
#' if(require("glmnet") && require("ranger")){
#' ## generate data
#' set.seed(1)
#' n  <- 200                                  # number of observations
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
#' learners <- c("lasso", "mlr3::lrn('ranger', num.trees = 30)")
#'
#' ## perform generic ML inference
#' # small number of splits to keep computation time low
#' x <- GenericML(Z, D, Y, learners, num_splits = 10,
#'                parallel = FALSE)
#'
#' ## access BLP generic targets for best learner w/o plot
#' get_BLP(x, learner = "best", plot = FALSE)
#'
#' ## access BLP generic targets for ranger learner w/o plot
#' get_BLP(x, learner = "mlr3::lrn('ranger', num.trees = 30)", plot = FALSE)
#'
#' ## access GATES generic targets for best learner w/o plot
#' get_GATES(x, learner = "best", plot = FALSE)
#'
#' ## access GATES generic targets for ranger learner w/o plot
#' get_GATES(x, learner = "mlr3::lrn('ranger', num.trees = 30)", plot = FALSE)
#'
#' ## access CLAN generic targets for "V1" & best learner, w/o plot
#' get_CLAN(x, learner = "best", variable = "V1", plot = FALSE)
#'
#' ## access CLAN generic targets for "V1" & ranger learner, w/o plot
#' get_CLAN(x, learner = "mlr3::lrn('ranger', num.trees = 30)",
#'          variable = "V1", plot = FALSE)
#' }
#'
#' @seealso
#' \code{\link{GenericML}},
#' \code{\link{get_BLP}},
#' \code{\link{get_CLAN}},
#' \code{\link[=plot.GenericML]{plot}}
#'
#' @export
get_GATES <- function(x, learner = "best", plot = TRUE){

  if(class(x) != "GenericML") stop("x needs to be instance of the class GenericML")

  if(learner == "best"){

    temp <- x$VEIN$best_learners$GATES

  } else{

    if(!(learner %in% x$arguments$learners_GenericML)){

      stop("Specified learner is not used in this instance of GenericML")

    } # IF

    temp <- x$VEIN$all_learners$GATES[[learner]]

  } # IF

  # only print the minimum of the two probabilities as p-value
  pval <- pmin(temp[, "Pr(<z) adjusted"], temp[, "Pr(>z) adjusted"])
  out  <- cbind(temp[, c("Estimate", "CB lower", "CB upper")],
                pval)
  colnames(out)[4] <- "Pr(>|z|)"

  if(plot) print(plot.GenericML(x = x, learner = learner, type = "GATES"))

  out

} # FUN


#' Accessor function for the CLAN generic target estimates
#'
#' @param x An object of the class \code{\link{GenericML}}.
#' @param variable The (character) name of a variabe on which CLAN was performed.
#' @param learner A character string of the learner whose CLAN generic target estimates shall be accessed. Default is \code{"best"} for the best learner for CLAN
#' @param plot Logical. If \code{TRUE} (default), a plot is printed.
#'
#' @return A numeric matrix of CLAN generic target estimates which contains information on point estimates, confidence bounds, and (adjusted) p-values. Furthermore, prints a plot if \code{plot = TRUE}.
#'
#' @examples
#' if(require("glmnet") && require("ranger")){
#' ## generate data
#' set.seed(1)
#' n  <- 200                                  # number of observations
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
#' learners <- c("lasso", "mlr3::lrn('ranger', num.trees = 30)")
#'
#' ## perform generic ML inference
#' # small number of splits to keep computation time low
#' x <- GenericML(Z, D, Y, learners, num_splits = 10,
#'                parallel = FALSE)
#'
#' ## access BLP generic targets for best learner w/o plot
#' get_BLP(x, learner = "best", plot = FALSE)
#'
#' ## access BLP generic targets for ranger learner w/o plot
#' get_BLP(x, learner = "mlr3::lrn('ranger', num.trees = 30)", plot = FALSE)
#'
#' ## access GATES generic targets for best learner w/o plot
#' get_GATES(x, learner = "best", plot = FALSE)
#'
#' ## access GATES generic targets for ranger learner w/o plot
#' get_GATES(x, learner = "mlr3::lrn('ranger', num.trees = 30)", plot = FALSE)
#'
#' ## access CLAN generic targets for "V1" & best learner, w/o plot
#' get_CLAN(x, learner = "best", variable = "V1", plot = FALSE)
#'
#' ## access CLAN generic targets for "V1" & ranger learner, w/o plot
#' get_CLAN(x, learner = "mlr3::lrn('ranger', num.trees = 30)",
#'          variable = "V1", plot = FALSE)
#' }
#'
#' @seealso
#' \code{\link{GenericML}},
#' \code{\link{get_BLP}},
#' \code{\link{get_GATES}},
#' \code{\link[=plot.GenericML]{plot}}
#'
#' @export
get_CLAN <- function(x, variable, learner = "best", plot = TRUE){

  if(class(x) != "GenericML") stop("x needs to be instance of the class GenericML")
  if(!(variable %in% names(x$VEIN$best_learners$CLAN))){
    stop(paste0("No CLAN was performed on this variable. ",
                "CLAN was performed on the variables ",
                paste0(names(x$VEIN$best_learners$CLAN), collapse = ", ")))
  }


  if(learner == "best"){

    temp <- x$VEIN$best_learners$CLAN[[variable]]

  } else{

    if(!(learner %in% x$arguments$learners_GenericML)){

      stop("Specified learner is not used in this instance of GenericML")

    } # IF

    temp <- x$VEIN$all_learners$CLAN[[learner]][[variable]]

  } # IF

  # only print the minimum of the two probabilities as p-value
  pval <- pmin(temp[, "Pr(<z) adjusted"], temp[, "Pr(>z) adjusted"])
  out  <- cbind(temp[, c("Estimate", "CB lower", "CB upper")],
                pval)
  colnames(out)[4] <- "Pr(>|z|)"

  if(plot) print(plot(x = x, learner = learner, type = "CLAN", CLAN_variable = variable))

  out

} # FUN

