#' Setup function for stratified sampling
#'
#' This function controls whether or not stratified sample splitting shall be performed. If no stratified sampling shall be performed, do not pass any arguments to this function (this is the default). If stratified sampling shall be performed, use this function to pass arguments to \code{\link[splitstackshape]{stratified}} in the package A character equal to covariance estimation function names in the \href{https://CRAN.R-project.org/package=splitstackshape}{splitstackshape}, which is used for obtaining the strata. In this case, the specification for \code{prop_aux} in \code{\link{GenericML}} does not have an effect because the number of samples in the auxiliary set is specified with the \code{size} argument in \code{\link[splitstackshape]{stratified}}.
#'
#' @param ... Named objects that shall be used as arguments in \code{\link[splitstackshape]{stratified}}. If empty (default), ordinary random sampling will be performed.
#'
#' @return
#' A list of named objects (possibly empty) specifying the stratified sampling strategy. If empty, no stratified sampling will be performed and instead ordinary random sampling will be performed.
#'
#' @details
#' The output of this setup function is intended to be used as argument \code{stratify} in the function \code{\link{GenericML}}. If arguments are passed to \code{\link[splitstackshape]{stratified}} via this function, make sure to  pass the necessary objects that \code{\link[splitstackshape]{stratified}} requires. The necessary objects are called \code{indt}, \code{group}, and \code{size} (see the documentation of  \code{\link[splitstackshape]{stratified}} for details). If either of these objects is missing, an error is thrown.
#'
#' @seealso
#' \code{\link[splitstackshape]{stratified}}
#' \code{\link{GenericML}}
#'
#' @examples
#' ## sample data of group membership (with two groups)
#' set.seed(1)
#' n <- 500
#' groups <- data.frame(group1 = rbinom(n, 1, 0.2),
#'                      group2 = rbinom(n, 1, 0.3))
#'
#' ## suppose we want both groups to be present in a strata...
#' group <- c("group1", "group2")
#'
#' ## ... and that the size of the strata equals half of the observations per group
#' size <- 0.5
#'
#' ## obtain a list of arguments that will be passed to splitstackshape::stratified()
#' setup_stratify(indt = groups, group = group, size = size)
#'
#' ## if no stratified sampling shall be used, do not pass anything
#' setup_stratify()
#'
#' @export
setup_stratify <- function(...) list(...)



#' Generates a function for sample splitting (internal use)
#'
#' Input checks will be performed via \code{\link{InputChecks_stratify}}.
#'
#' @param args_stratified A list of arguments that shall be passed to \code{\link[splitstackshape]{stratified}}; typically returned by \code{\link{setup_stratify}}
#' @param N Number of samples
#' @param prop_aux Proportion of samples that shall be in the auxiliary set. In case of stratified sampling, the effective proportion of samples in the auxiliary set might deviate from \code{prop_aux}, depending on the specifications of the strata.
#'
#' @return
#' A function without arguments that performs random sample splitting, either stratified or non-stratified. Specifically, it returns indices of samples in the auxiliary set.
#'
#' @noRd
make_split_fn <- function(args_stratified, N, prop_aux)
{

  if(length(args_stratified) == 0L)
  {
    ## case 1: no stratified sampling shall be used, so sample randomly
    fn <- function() sample(x = 1:N, size = floor(prop_aux * N), replace = FALSE)
    fn <- structure(fn, type = "nostratification")

  } else
  {
    ## case 2: stratified sampling
    # overwrite/specify the arguments 'keep.rownames' and 'bothSets'
    args_stratified$keep.rownames <- TRUE
    args_stratified$bothSets      <- FALSE

    # make 'indt' a data frame with row names because we want to retain
    # row indices, and data tables don't support row indices
    indt            <- as.data.frame(args_stratified$indt)
    if(nrow(indt) != N) stop("'indt' must have the same number of rows as Z")
    rownames(indt)  <- 1:N
    args_stratified$indt <- indt

    # return the function
    fn <- function()
    {
      out <- do.call(what = splitstackshape::stratified,
                     args = args_stratified)
      as.integer(out$rn)
    } # FUN

    fn <- structure(fn, type = "stratified")

  } # IF

  return(fn)

} # FUN


#' Performs the sample splitting (internal use)
#'
#' @param split_fn A function that splits the data by returning the indices of samples that will be in the auxiliary set. Such a function is returned by \code{make_split_fn} (internal function).
#' @param D Binary vector of treatment assignment
#' @param N sample size, must be equal to length of \code{D}
#'
#' @return
#' A list consisting of the indices of samples in the main set (\code{M_set}) and auxiliary set (\code{A_set}).
#'
#' @noRd
sample_split <- function(split_fn, D, N){

  # initialize
  temp <- TRUE

  while(temp){

    # sample candidate set for A_set
    A_set <- split_fn()

    # Avoid imbalance in A_set for estimation of BCA and CATE.
    # BCA is estimated on the control units, CATE on the treated units.
    # We want to avoid that either of them is estimated on too small a sample.
    # To achieve this, have the control units in A_set make up no more
    # than 90% of all samples in A_set
    if(mean(D[A_set] == 0) <= 0.9) temp <- FALSE

  } # WHILE

  A_set <- as.integer(sort(A_set, decreasing = FALSE))

  # return
  return(list(A_set = A_set,
              M_set = setdiff(1:N, A_set)))

} # FUN
