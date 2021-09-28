print.GenericML <- function(GenericML.obj){
  
  cat("GenericML object with the following specifications:\n")
  cat("\t* Propensity Score learner:", GenericML.obj$arguments$learner.propensity.score, "\n")
  cat("\t* Generic ML learners:", paste(GenericML.obj$arguments$learners.genericML, collapse = ", "), "\n")
  cat("\t* S =", GenericML.obj$arguments$num.splits, "splits are used\n")
  cat("\t*", ifelse(GenericML.obj$arguments$HT.transformation, "A", "No"), "HT transformation is used\n")
  cat("\nThe", 100 * (1-2*GenericML.obj$arguments$significance.level), "% confidence bounds of the best BLP estimates are given by\n")
  cat("\t beta.1: (", 
      round(GenericML.obj$VEIN$best.learners$BLP["beta.1", "CB lower"], 3), ",",
      round(GenericML.obj$VEIN$best.learners$BLP["beta.1", "CB upper"], 3), ")")
  cat("\t beta.2: (", 
      round(GenericML.obj$VEIN$best.learners$BLP["beta.2", "CB lower"], 3), ",",
      round(GenericML.obj$VEIN$best.learners$BLP["beta.2", "CB upper"], 3), ")\n")   
  cat("The best learner for the CATE is ", GenericML.obj$best.learners$best.learner.for.CATE,
      " (lambda of ", round(max(GenericML.obj$best.learners$lambda.overview[,"lambda"]), 3), ")\n",
      "The best learner for the GATES and CLAN is ", GenericML.obj$best.learners$best.learner.for.GATES,
      " (lambda.bar of ", round(max(GenericML.obj$best.learners$lambda.overview[,"lambda.bar"]), 3), ")", sep = "")
  
} # FUN
