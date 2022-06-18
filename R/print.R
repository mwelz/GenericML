#' Print method for a \code{GenericML} object
#'
#' Prints key results of the analyses conducted in \code{\link{GenericML}()}.
#'
#' @param x An object of the class \code{"\link{GenericML}"}, as returned by the function \code{\link{GenericML}()}.
#' @param digits Number of digits to print.
#' @param ... Additional arguments to be passed down.
#'
#' @return
#' A print to the console.
#'
#' @examples
#' if(require("ranger")){
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
#' ## specify learners
#' learners <- c("random_forest")
#'
#' ## perform generic ML inference
#' # small number of splits to keep computation time low
#' x <- GenericML(Z, D, Y, learners, num_splits = 2,
#'                parallel = FALSE)
#'
#' ## print
#' print(x)
#' }
#'
#' @export
print.GenericML <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  if(!inherits(x = x, what = "GenericML", which = FALSE)){
    stop("The object 'x' must be an instance of GenericML()")
  }

  if(is.numeric(x$arguments$learner_propensity_score)){
    prop.lrn <- "user-supplied"
  } else{
    prop.lrn <- x$arguments$learner_propensity_score
  }

  cat("GenericML object with the following specifications:\n")
  cat("\t* Propensity score learner:", prop.lrn, "\n")
  cat("\t* Generic ML learners:", paste(x$arguments$learners_GenericML, collapse = ", "), "\n")
  cat("\t* S =", x$arguments$num_splits, "splits are used\n")
  cat("\t*", ifelse(x$arguments$HT, "A", "No"), "HT transformation is used\n")
  cat("\nThe", 100 * (1-2*x$arguments$significance_level), "% confidence intervals of the best BLP estimates are given by\n")
  cat("\t beta.1: (",
      round(x$VEIN$best_learners$BLP["beta.1", "CB lower"], 3), ",",
      round(x$VEIN$best_learners$BLP["beta.1", "CB upper"], 3), ")")
  cat("\t beta.2: (",
      round(x$VEIN$best_learners$BLP["beta.2", "CB lower"], 3), ",",
      round(x$VEIN$best_learners$BLP["beta.2", "CB upper"], 3), ")\n")
  cat("The best learner for the BLP is ", x$best$BLP,
      " (lambda of ", round(max(x$best$overview[,"lambda"]), digits), ")\n",
      "The best learner for the GATES and CLAN is ", x$best$GATES,
      " (lambda.bar of ", round(max(x$best$overview[,"lambda.bar"]), digits), ")", sep = "")
  cat("\n")

} # FUN


#' Print method for a \code{"BLP_info"} object
#'
#' @param x An object of the class \code{"BLP_info"}, as returned by the function \code{\link{get_BLP}()}.
#' @param digits Number of digits to print.
#' @param ... Additional arguments to be passed down.
#'
#' @return
#' A print to the console.
#'
#' @export
print.BLP_info <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  if(!inherits(x = x, what = "BLP_info", which = FALSE)){
    stop("The object 'x' must be of class BLP_info")
  }

  cat("BLP generic targets\n---\n")
  print_accessors(x = x, digits = digits, print_confidence = TRUE, ... = ...)
} # FUN


#' Print method for a \code{"GATES_info"} object
#'
#' @param x An object of the class \code{"GATES_info"}, as returned by the function \code{\link{get_GATES}()}.
#' @param digits Number of digits to print.
#' @param ... Additional arguments to be passed down.
#'
#' @return
#' A print to the console.
#'
#' @export
print.GATES_info <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  if(!inherits(x = x, what = "GATES_info", which = FALSE)){
    stop("The object 'x' must be of class GATES_info")
  }

  cat("GATES generic targets\n---\n")
  print_accessors(x = x, digits = digits, print_confidence = TRUE, ... = ...)
} # FUN


#' Print method for a \code{"CLAN_info"} object
#'
#' @param x An object of the class \code{"CLAN_info"}, as returned by the function \code{\link{get_CLAN}()}.
#' @param digits Number of digits to print.
#' @param ... Additional arguments to be passed down.
#'
#' @return
#' A print to the console.
#'
#' @export
print.CLAN_info <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  if(!inherits(x = x, what = "CLAN_info", which = FALSE)){
    stop("The object 'x' must be of class CLAN_info")
  }

  cat(sprintf("CLAN generic targets for variable '%s' \n---\n", x$CLAN_variable))
  print_accessors(x = x, digits = digits, print_confidence = TRUE, ... = ...)
} # FUN


#' Internal function to support print methods of accessors; returns a print
#'
#' @param x An object of class \code{"BLP_info"}, \code{"GATES_info"}, or \code{"CLAN_info"}
#' @param digits Number of digits to print
#' @param print_confidence Logical. Shall the confidence level be printed as well?
#' @param ... Additional arguments to be passed down
#'
#' @noRd
print_accessors <- function(x, digits, print_confidence,...)
{
  # prepare coefficient matrix
  mat <- cbind(x$estimate, x$confidence_interval, x$p_value)
  colnames(mat) <- c("Estimate", "CI lower", "CI upper", "p value")

  # print it
  stats::printCoefmat(mat, digits = digits, ... = ...)

  # print confidence level if requested
  if(isTRUE(print_confidence))
  {
    cat("---\nConfidence level of confidence interval [CI lower, CI upper]: ",
        format(100 * x$confidence_level), " %\n", sep = "")
  } # IF

  # create plot if requested
  p <- x$plot
  if (!is.null(p)) print(p)

  # return object invisibly
  invisible(x)
} # FUN


#' Print method for a \code{"heterogeneity_CLAN"} object
#'
#' @param x An object of class \code{"\link{heterogeneity_CLAN}"}.
#' @param ... Additional arguments to be passed down.
#'
#' @return
#' A print to the console.
#'
#' @export
print.heterogeneity_CLAN <- function(x, ...)
{
  if(!inherits(x = x, what = "heterogeneity_CLAN", which = FALSE)){
    stop("x needs to be instance of the class heterogeneity_CLAN")
  } # IF

  cat("There is significant treatment effect heterogeneity along the following CLAN variables:\n\t")

  if(x$significant$num_variables > 0L)
  {
    cat(x$significant$variables, sep = ", ")
  } else{
    cat("(none)")
  } # IF

  cat("\n---\nLevel of significance: ",
      format(100 * x$significance_level), " %\n", sep = "")
  cat("---\nUse 'get_CLAN()' to further explore CLAN along these variables\n")

  # return object invisibly
  invisible(x)
} # FUN


#' @export
print.best <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  print(x$overview, digits = digits)
  cat("---\n")
  cat("The best learner for BLP is ", x$BLP,
      " with lambda = ", round(x$overview[x$BLP, "lambda"], digits), ".\n", sep = "")
  cat("The best learner for GATES and CLAN is ", x$BLP,
      " with lambda.bar = ", round(x$overview[x$GATES, "lambda.bar"], digits),
      ".\n", sep = "")
}
