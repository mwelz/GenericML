#' Setup function for stratified sampling
#'
#' Stratified sampling will be performed via \code{\link[splitstackshape]{stratified}}. This function gathers the arguments that shall be passed to \code{\link[splitstackshape]{stratified}}. If no arguments are passed, no stratified sampling will be performed. Instead, ordinary random sampling will be performed.
#'
#' @param ... Named objects that shall be used as arguments in \code{\link[splitstackshape]{stratified}}. If empty, ordinary random sampling will be performed.
#'
#' @return
#' A list of named objects (possibly empty) specifying the stratified sampling strategy. If empty, no stratified sampling will be performed. Instead, ordinary random sampling will be performed.
#'
#' @export
setup_stratified <- function(...) list(...)



#' Generates a function for sample splitting (internal use)
#'
#' @param args_stratified A list of arguments that shall be passed to \code{\link[splitstackshape]{stratified}}; typically returned by \code{\link{setup_stratified}}
#' @param N Number of samples
#' @param prop_aux Proportion of samples that shall be in the auxiliary set. In case of stratified sampling, the effective proportion of samples in the auxiliary set might deviate from \code{prop_aux}, depending on the specifications of the strata.
#'
#' @return
#' A function without arguments that performs random sample splitting, either stratified or non-stratified. Specifically, it returns indices of samples in the auxiliary set.
#'
#' @noRd
make_split_fn <- function(args_stratified, N, prop_aux)
{

  ## ensure that 'args_stratified' is a list
  if(!is.list(args_stratified)){
    stop("'args_stratified' must be a list, for instance as returned by setup_stratified()")
  }

  if(length(args_stratified) == 0L)
  {
    ## case 1: no stratified sampling shall be used, so sample randomly
    fn <- function() sample(x = 1:N, size = floor(prop_aux * N), replace = FALSE)
    fn <- structure(fn, type = "nostratification")

  } else
  {
    ## case 2: stratified sampling
    # check that all necessary arguments for splitstackshape::stratified are passed
    if(!all(c("indt", "group", "size") %in% names(args_stratified))){
      stop(paste0("splitstackshape::stratified requires at least the arguments ",
                  "'indt', 'group', and 'size', which were not passed to setup_stratified().",
                  " See ?splitstackshape::stratified for details." ))
    } # IF

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
#' @param D Binary vector of treatment assignment
#' @param N sample size
#' @param N_set 1:N
#' @param prop Total number of samples in the auxiliary set
#'
#' @noRd
sample_split <- function(D, N, N_set = 1:N, prop){

  temp <- TRUE

  while(temp){

    # sample candidate set for A_set
    A_set <- sample(N_set, size = prop, replace = FALSE)

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
              M_set = setdiff(N_set, A_set)))

} # FUN
