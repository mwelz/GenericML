#' Print method for a \code{GenericML} object
#'
#' Prints key results of the analyses conducted in \code{\link{GenericML}}.
#'
#' @param x An instance of \code{\link{GenericML}}.
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
print.GenericML <- function(x, ...){

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
  cat("\nThe", 100 * (1-2*x$arguments$significance_level), "% confidence bounds of the best BLP estimates are given by\n")
  cat("\t beta.1: (",
      round(x$VEIN$best_learners$BLP["beta.1", "CB lower"], 3), ",",
      round(x$VEIN$best_learners$BLP["beta.1", "CB upper"], 3), ")")
  cat("\t beta.2: (",
      round(x$VEIN$best_learners$BLP["beta.2", "CB lower"], 3), ",",
      round(x$VEIN$best_learners$BLP["beta.2", "CB upper"], 3), ")\n")
  cat("The best learner for the BLP is ", x$best$BLP,
      " (lambda of ", round(max(x$best$overview[,"lambda"]), 3), ")\n",
      "The best learner for the GATES and CLAN is ", x$best$GATES,
      " (lambda.bar of ", round(max(x$best$overview[,"lambda.bar"]), 3), ")", sep = "")
  cat("\n")

} # FUN
