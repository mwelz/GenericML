#' prints information on an instance of \code{\link{GenericML}}.
#'
#' @param x An instance of \code{\link{GenericML}}.
#' @param ... additional arguments to be passed
#'
#' @export
print.GenericML <- function(x, ...){

  cat("GenericML object with the following specifications:\n")
  cat("\t* Propensity score learner:", x$arguments$learner_propensity_score, "\n")
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

} # FUN
