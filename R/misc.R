
#' Find lower and upper median of a vector.
#'
#' See Comment 4.2 in the paper.
#'
#' @param x A numeric vector whose medians we are interested in.
#'
#' @return A list with the upper and lower median and the Med statistic.
#'
#' @export
Med <- function(x){

  # get the empirical CDF of X
  ecdf.x <- stats::ecdf(x)

  # evaluate the ensuing probabilities
  F.x <- ecdf.x(x)

  # get lower median and upper median
  lower       <- min(x[F.x >= 0.5])
  upper.array <- x[(1 - F.x) >= 0.5]
  upper       <- ifelse(length(upper.array) == 0,
                        lower,
                        max(upper.array)) # account for case where upper.array is empty

  return(list(lower_median = lower,
              upper_median = upper,
              Med = mean(c(lower, upper))))

} # END FUN



#' Partition a vector into groups based on its quantiles.
#'
#' Partitions a vector into quantile groups and returns a logical matrix indicating group membership.
#'
#' @param x The vector to be partitioned
#' @param cutoffs The quantile cutoffs for the partition. Default are the quartiles: \code{c(0.25, 0.5, 0.75)}.
#' @param names_quantile Logical. If \code{TRUE}, then the column names of the returned matrix are the quantiles as in \code{cutoffs}. If \code{FALSE}, the names are the numeric intervals that constitute the grouping.
#'
#' @return An object of the class \code{quantile_group}, which is a logical matrix indicating group membership
#'
#' @export
quantile_group <- function(x,
                           cutoffs = c(0.25, 0.5, 0.75),
                           names_quantile = TRUE){
  # cutoffs are the quantile cutoffs (like c(0.25, 0.5, 0.75))

  # get quatiles
  q         <- stats::quantile(x, cutoffs)
  q         <- c(-Inf, q, Inf)

  # check if breaks are unique: if x exhibits low variation, there might be empty quantile bins, which can cause an error in the cut() function. In this case, we add random noise to x to induce variation. NB: this bug has been spotted and fixed by Lucas Kitzmueller. All credits for this fix go to him!
  if(length(unique(q)) != length(q)){

    # specify standard deviation of the noise (x may have zero variation)
    sd <- ifelse(stats::var(x) == 0, 0.001, sqrt(stats::var(x) / 20))

    # add noise and updare quantiles
    x <- x + stats::rnorm(length(x), mean = 0, sd = sd)
    q <- stats::quantile(x, cutoffs)
    q <- c(-Inf, q, Inf)

  } # IF

  groups    <- as.character(cut(x, breaks = q, include.lowest = TRUE, right = FALSE, dig.lab = 3))
  group.nam <- unique(groups)
  group.nam <- group.nam[order(
    as.numeric(substr(sub("\\,.*", "", group.nam), 2, stop = 1e8L)),
    decreasing = FALSE)] # ensure the order is correct
  group.mat <- matrix(NA, length(x), length(group.nam))
  nam       <- rep(NA, length(group.nam))

  for(j in 1:length(group.nam)){
    if(j == 1){
      nam[j] <- paste0("<", 100*cutoffs[j], "% quantile")
    } else if (j == length(group.nam)){
      nam[j] <- paste0(">=", 100*cutoffs[j-1], "% quantile")
    } else{
      nam[j] <- paste0("[", 100*cutoffs[j-1], ",", 100*cutoffs[j], ")% quantile")
    }
    group.mat[,j] <- groups == group.nam[j]
  }

  if(names_quantile){
    colnames(group.mat) <- nam
  } else{
    colnames(group.mat) <- group.nam
  }
  return(structure(group.mat, type = "quantile_group"))
} # FUN



#' Performs generic ML for a given learning technique and a given split of the data.
#'
#' Can be seen as a single iteration of Algorithm 1 in the paper.
#'
#' @param Z A matrix or data frame of the covariates.
#' @param D A binary vector of treatment assignment.
#' @param Y The response vector.
#' @param learner A string specifying the machine learner to be used for estimating the BCA and CATE. Either \code{'elastic_net'}, \code{'random_forest'}, or \code{'tree'}. Can alternatively be specified by using \code{mlr3} syntax \emph{without} specification if the learner is a regression learner or classification learner. Example: \code{'mlr3::lrn("ranger", num.trees = 500)'} for a random forest learner. Note that this is a string and the absence of the \code{classif.} or \code{regr.} keywords. See \url{https://mlr3learners.mlr-org.com} for a list of \code{mlr3} learners.
#' @param propensity_scores A numeric vector of propensity scores.
#' @param M_set a numerical vector of indices of observations in the main sample.
#' @param A_set a numerical vector of indices of observations in the auxiliary sample. Default is complementary set to \code{M_set}.
#' @param Z_CLAN A matrix of variables that shall be considered for the CLAN. Each column represents a variable for which CLAN shall be performed. If \code{NULL} (default), then \code{Z_CLAN = Z}, i.e. CLAN is performed for all variables in \code{Z}.
#' @param HT Logical. If \code{TRUE}, a HT transformation is applied in the BLP and GATES regressions. Default is \code{FALSE}.
#' @param quantile_cutoffs The cutoff points of quantiles that shall be used for GATES grouping. Default is \code{c(0.25, 0.5, 0.75)}, which corresponds to the four quartiles.
#' @param X1_BLP Specifies the design matrix \eqn{X_1} in the BLP regression. See the documentation of \code{\link{setup_X1}} for details.
#' @param X1_GATES Same as \code{X1_BLP}, just for the the GATES regression.
#' @param diff_GATES Specifies the generic targets of GATES. See the documentation of \code{\link{setup_diff}} for details.
#' @param diff_CLAN Same as \code{diff_GATES}, just for the CLAN generic targets.
#' @param vcov_BLP Specifies the covariance matrix estimator in the BLP regression. See the documentation of \code{\link{setup_vcov}} for details.
#' @param vcov_GATES Same as \code{vcov_BLP}, just for the GATES regression.
#' @param equal_variances_CLAN Logical. If \code{TRUE}, the the two within-group variances of the differences between the CLAN generic targets are assumed to be equal. Default is \code{FALSE}.
#' @param significance_level Significance level for VEIN. Default is 0.05.
#' @param min_variation Minimum variation of the predictions before random noise with distribution \eqn{N(0, var(Y)/20)} is added. Default is \code{1e-05}.
#'
#' @return a list with instances of the classes \code{BLP}, \code{GATES}, \code{CLAN}, \code{proxy_BCA}, and \code{proxy_CATE}. In addition, the lambda parameters for finding the best learner are returned.
#'
#' @export
GenericML_single <- function(Z, D, Y,
                             learner,
                             propensity_scores,
                             M_set,
                             A_set                = setdiff(1:length(Y), M_set),
                             Z_CLAN               = NULL,
                             HT                   = FALSE,
                             quantile_cutoffs     = c(0.25, 0.5, 0.75),
                             X1_BLP               = setup_X1(),
                             X1_GATES             = setup_X1(),
                             diff_GATES           = setup_diff(),
                             diff_CLAN            = setup_diff(),
                             vcov_BLP             = setup_vcov(),
                             vcov_GATES           = setup_vcov(),
                             equal_variances_CLAN = FALSE,
                             significance_level   = 0.05,
                             min_variation        = 1e-05){

  # input checks
  InputChecks_D(D)
  InputChecks_Y(Y)
  InputChecks_Z(Z)
  InputChecks_equal.length3(D, Y, Z)
  InputChecks_X1(X1_BLP)
  InputChecks_X1(X1_GATES)
  InputChecks_vcov.control(vcov_BLP)
  InputChecks_vcov.control(vcov_GATES)
  InputChecks_diff(diff_GATES, K = length(quantile_cutoffs) + 1)
  InputChecks_diff(diff_CLAN, K = length(quantile_cutoffs) + 1)

  if(is.null(Z_CLAN)) Z_CLAN <- Z # if no input provided, set it equal to Z
  learner <- get.learner_regr(make.mlr3.environment(learner, regr = TRUE))


  # call the main function
  GenericML_single_NoChecks(Z = Z, D = D, Y = Y,
                            propensity_scores = propensity_scores,
                            learner = learner,
                            M_set = M_set, A_set = A_set,
                            Z_CLAN                     = Z_CLAN,
                            X1_BLP                     = X1_BLP,
                            X1_GATES                   = X1_GATES,
                            HT                         = HT,
                            vcov_BLP                   = vcov_BLP,
                            vcov_GATES                 = vcov_GATES,
                            equal_variances_CLAN       = equal_variances_CLAN,
                            quantile_cutoffs           = quantile_cutoffs,
                            diff_GATES                 = diff_GATES,
                            diff_CLAN                  = diff_CLAN,
                            significance_level         = significance_level,
                            min_variation              = min_variation)

} # END FUN


# helper that skips the input checks
GenericML_single_NoChecks <-
  function(Z, D, Y,
           propensity_scores,
           learner,
           M_set, A_set,
           Z_CLAN                     = NULL,
           X1_BLP                     = setup_X1(),
           X1_GATES                   = setup_X1(),
           HT                         = FALSE,
           vcov_BLP                   = setup_vcov(),
           vcov_GATES                 = setup_vcov(),
           equal_variances_CLAN       = FALSE,
           quantile_cutoffs           = c(0.25, 0.5, 0.75),
           diff_GATES                 = setup_diff(),
           diff_CLAN                  = setup_diff(),
           significance_level         = 0.05,
           min_variation              = 1e-05){


    ### step 1a: learn proxy predictors by using the auxiliary set ----

    # get proxy baseline estimates
    proxy_BCA.obj <- proxy_BCA_NoChecks(
      Z = Z, D = D, Y = Y,
      A_set         = A_set,
      learner       = learner,
      min_variation = min_variation)

    # get estimates on main sample
    proxy_BCA_M     <- proxy_BCA.obj$estimates[M_set]

    # get the proxy estimates of the CATE
    proxy_CATE.obj <-
      proxy_CATE_NoChecks(
        Z = Z, D = D, Y = Y,
        A_set          = A_set,
        learner        = learner,
        proxy_BCA = proxy_BCA.obj$estimates,
        min_variation  = min_variation)

    # get estimates on main sample
    proxy_CATE_M <- proxy_CATE.obj$estimates$CATE[M_set]


    ### step 1b: estimate BLP parameters ----
    blp.obj <- BLP_NoChecks(
      D                  = D[M_set],
      Y                  = Y[M_set],
      propensity_scores  = propensity_scores[M_set],
      proxy_BCA     = proxy_BCA_M,
      proxy_CATE         = proxy_CATE_M,
      HT                 = HT,
      X1_control         = setup_X1(),
      vcov_control       = vcov_BLP,
      significance_level = significance_level)



    ### step 1c: estimate GATES parameters by OLS ----
    # group the proxy estimators for the CATE in the main sample by quantiles
    membership_M <- quantile_group(proxy_CATE_M,
                                   cutoffs = quantile_cutoffs,
                                   names_quantile = TRUE)

    # estimate GATES
    gates.obj <- GATES_NoChecks(
      D                   = D[M_set],
      Y                   = Y[M_set],
      propensity_scores   = propensity_scores[M_set],
      proxy_BCA      = proxy_BCA_M,
      proxy_CATE          = proxy_CATE_M,
      membership          = membership_M,
      HT                  = HT,
      X1_control          = setup_X1(),
      vcov_control        = vcov_GATES,
      diff                = diff_GATES,
      significance_level  = significance_level)


    ### step 1d: estimate CLAN parameters on the main sample ----
    clan.obj <- CLAN_NoChecks(
      Z_CLAN             = Z_CLAN[M_set,,drop = FALSE],
      membership         = membership_M,
      equal_variances    = equal_variances_CLAN,
      diff               = diff_CLAN,
      significance_level = significance_level)


    ### step 1e: get parameters over which we maximize to find the "best" ML method ----
    best.obj <- lambda_parameters(BLP        = blp.obj,
                                  GATES      = gates.obj,
                                  proxy_CATE = proxy_CATE_M,
                                  membership = membership_M)

    ### organize output in a list ----
    return(list(BLP            = blp.obj,
                GATES          = gates.obj,
                CLAN           = clan.obj,
                proxy_CATE     = proxy_CATE.obj,
                proxy_BCA = proxy_BCA.obj,
                best           = best.obj
    ))

  } # FUN



#' Returns the two lambda parameters that are used to find the best ML method.
#'
#' @param BLP An instance of the class \code{BLP}.
#' @param GATES An instance of the class \code{GATES}.
#' @param proxy_CATE Proxy estimates of the CATE.
#' @param membership A logical matrix that indicates the group membership of each observation. Needs to be an instance of \code{\link{quantile_group}}.
#'
#' @return A list containing the parameters \code{lambda} and \code{lambda.bar}.
#'
#' @export
lambda_parameters <- function(BLP,
                              GATES,
                              proxy_CATE,
                              membership){

  temp <- GATES$coefficients
  gates.coefs <- temp[startsWith(rownames(temp), "gamma."), "Estimate"]

  return(list(lambda = BLP$coefficients["beta.2", "Estimate"]^2 * stats::var(proxy_CATE),
              lambda.bar = as.numeric(colMeans(membership) %*%  gates.coefs^2)))

} # END FUN


#' performs stratified sampling
#' @param group The column or columns that should be used to create the groups. Can be a character vector of column names (recommended) or a numeric vector of column.
#' @param relative.size proportional to the number of observations per group.
#' @param select A named list containing levels from the "group" variables in which you are interested. The list names must be present as variable names for the input dataset.
#'
#' #################################################
#' The function below effectively does the following:
#' DF <- data.frame(
# ID = 1:100,
# A = sample(c("AA", "BB", "CC", "DD", "EE"), 100, replace = TRUE),
# B = rnorm(100), C = abs(round(rnorm(100), digits=1)),
# D = sample(c("CA", "NY", "TX"), 100, replace = TRUE),
# E = sample(c("M", "F"), 100, replace = TRUE))
#
# group = DF[,c("E", "D")] # user supplied. This matrix/vector needs to contain group values
# select = list(E = "M", D = c("CA", "TX")) # user supplied with levels
#
# group <- as.matrix(group) # TODO: if select doesn't specify levels?
#
# eligible.samples <- lapply(colnames(group), function(j) which(group[,j] %in% select[[j]]) ) # select samples that satisfy the selection criteria within each group
# pool <- sort(Reduce(intersect, eligible.samples)) # find intersection of groupwise eligible samples
# sample(pool, floor(length(pool) * relative.size)) # sample from the pool
#
#  !!!!!!!!! there is a number of potential issues here: !!!!!!!!!
# 1. it holds that |stratified samples| <= nrow(group) * relative size. Thus, the strata can become VERY small, which is not optimal for fitting a ML method.
# 2. statistical concern: this is nonrandom sample selection, hence it might induce a sampling bias.
# 3. What shall the role of the strata be? Shall it be the A set or the M set? Suppose it is the A set, is the A set then [n]\A or pool\A? The second option implies that we effectively only work on the pool.
#'
#' @noRd
stratified <- function(group, relative.size, select = NULL){

  # perform stratified sampling
  smpl <- splitstackshape::stratified(indt = group, group = group, size = relative.size, select = select)
  return(sort(smpl$ID))

} # FUN
