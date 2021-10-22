#' Accessor function for the BLP generic target estimates after VEIN
#'
#' @param x An object of the class \code{\link{GenericML}}.
#' @param learner A character string of the learner whose BLP generic target estimates shall be accessed. Default is \code{"best"} for the best learner for BLP.
#' @param plot Logical. If \code{TRUE} (default), a plot is printed.
#'
#' @return A numeric matrix of BLP generic target estimates
#'
#' @export
get_BLP <- function(x, learner = "best", plot = TRUE){

  if(class(x) != "GenericML") stop("x needs to be instance of the class GenericML")

  if(learner == "best"){

    out <- x$VEIN$best_learners$BLP

  } else{

    if(!(learner %in% x$arguments$learners_GenericML)){

      stop("Specified learner is not used in this instance of GenericML")

    } # IF

    out <- x$VEIN$all_learners$BLP[[learner]]

  } # IF

  if(plot) print(plot.GenericML(x = x, learner = learner, type = "BLP"))

  out

} # FUN




#' Accessor function for the GATES generic target estimates after VEIN
#'
#' @param x An object of the class \code{\link{GenericML}}.
#' @param learner A character string of the learner whose GATES generic target estimates shall be accessed. Default is \code{"best"} for the best learner for GATES.
#' @param plot Logical. If \code{TRUE} (default), a plot is printed.
#'
#' @return A numeric matrix of GATES generic target estimates
#'
#' @export
get_GATES <- function(x, learner = "best", plot = TRUE){

  if(class(x) != "GenericML") stop("x needs to be instance of the class GenericML")

  if(learner == "best"){

    out <- x$VEIN$best_learners$GATES

  } else{

    if(!(learner %in% x$arguments$learners_GenericML)){

      stop("Specified learner is not used in this instance of GenericML")

    } # IF

    out <- x$VEIN$all_learners$GATES[[learner]]

  } # IF

  if(plot) print(plot.GenericML(x = x, learner = learner, type = "GATES"))

  out

} # FUN


#' Accessor function for the CLAN generic target estimates after VEIN
#'
#' @param x An object of the class \code{\link{GenericML}}.
#' @param variable The (character) name of a variabe on which CLAN was performed.
#' @param learner A character string of the learner whose CLAN generic target estimates shall be accessed. Default is \code{"best"} for the best learner for CLAN
#' @param plot Logical. If \code{TRUE} (default), a plot is printed.
#'
#' @return A numeric matrix of CLAN generic target estimates
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

    out <- x$VEIN$best_learners$CLAN[[variable]]

  } else{

    if(!(learner %in% x$arguments$learners_GenericML)){

      stop("Specified learner is not used in this instance of GenericML")

    } # IF

    out <- x$VEIN$all_learners$CLAN[[learner]][[variable]]

  } # IF

  if(plot) print(plot(x = x, learner = learner, type = "CLAN", CLAN_variable = variable))

  out

} # FUN

