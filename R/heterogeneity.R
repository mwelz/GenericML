#' Evaluate treatment effect heterogeneity along CLAN variables
#'
#' This function tests for statistical significance of all CLAN difference parameters that were specified in the function \code{\link{setup_diff}()}. It reports all CLAN variables along which there are significant difference parameters, which corresponds to evidence for treatment effect heterogeneity along this variable, at the specified significance level.
#'
#' @param x An object of class \code{"\link{GenericML}"}, as returned by the function \code{\link{GenericML}()}.
#' @param learner A character string of the learner whose CLAN generic target estimates are of interest. Default is \code{"best"} for the best learner for CLAN.
#' @param significance_level Level for the significance tests. Default is 0.05.
#'
#' @return
#' An object of class \code{"heterogeneity_CLAN"}, consisting of the following components:
#' \describe{
#'   \item{\code{p_values}}{A matrix of p values of all CLAN difference parameters for all CLAN variables.}
#'   \item{\code{significant}}{The names of variables with at least one significant CLAN difference parameter (\code{"variables"}), their number \code{"num_variables"}, and the total number of significant CLAN difference parameters \code{"num_params"}. All significance tests were performed at level \code{significance_level}.}
#'   \item{\code{min_pval}}{Information on the smallest p value: Its value (\code{"value"}), the variable in which it was estimated (\code{"variable"}), the CLAN difference parameter it belongs to (\code{"parameter"}), and whether or not it is significant at level \code{significance_level} (\code{"significant"}).}
#'   \item{\code{"learner"}}{Name of the learner whose median estimates we used for the listed results.}
#'   \item{\code{"significance_level"}}{The level of the significance tests.}
#'   }
#'
#' @seealso
#' \code{\link{GenericML}()}
#'
#' @export
heterogeneity_CLAN <- function(x, learner = "best", significance_level = 0.05)
{
  # input check
  isGenericMLcheck(x)
  stopifnot(is.numeric(significance_level) && length(significance_level) == 1L)
  stopifnot(0.0 < significance_level && significance_level < 0.5)

  # get CLAN names
  variables <- CLAN_names(x)

  # get number of groups
  K   <- length(x$arguments$quantile_cutoffs) + 1L
  toK <- seq_len(K)

  # get CLAN parameters
  CLAN_params <- rownames(x$VEIN$best_learners$CLAN[[1L]])

  # prepare matrix of p values
  pvalues_mat <- matrix(NA_real_, nrow = length(CLAN_params) - K,
                        ncol = length(variables))
  rownames(pvalues_mat) <- CLAN_params[-toK]
  colnames(pvalues_mat) <- variables


  for(variable in variables)
  {
    # access the results
    results <- accessor_CLAN_noChecks(x = x,
                                      variable = variable,
                                      learner = learner)

    # get the adjusted p values
    pvalues <- accessor_output(x = x, accessor_obj = results, plot = FALSE,
                               type = "CLAN", learner = learner,
                               CLAN_variable = variable, ATE = FALSE)$p_value

    # get p values of the difference parameters
    pvalues_mat[,variable] <- pvalues[-toK]

  } # FOR

  # evaluate significance
  pvalues_mat_sig <- pvalues_mat < significance_level

  # variables which had at least one significant difference parameter
  sig_vars <- variables[colSums(pvalues_mat_sig) > 0L]

  # number of significant difference parameters
  num_sig_params <- sum(rowSums(pvalues_mat_sig))

  # smallest p value of difference parameters
  minpval <- min(pvalues_mat)
  idx <-
    which(pvalues_mat == minpval, arr.ind = TRUE)[1L,]
  minpval_param <- rownames(pvalues_mat)[idx["row"]] # difference parameter with minimum p value
  minpval_var   <- variables[idx["col"]] # variable with minimum p value
  significant   <- minpval < significance_level # is it significant?

  return(
    structure(list(p_values = pvalues_mat,
                   significant = list(variables = sig_vars,
                                      num_variables = length(sig_vars),
                                      num_params = num_sig_params),
                   min_pval = list(value = minpval, variable = minpval_var,
                                   parameter = minpval_param, significant = significant),
                   learner = learner, significance_level = significance_level),
              class = "heterogeneity_CLAN")
  )

} # FUN
