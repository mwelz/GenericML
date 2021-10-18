
# get lower and upper median as in comment 4.2 in Chernozhukov et al. (2021)
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

  return(list(lower.median = lower,
              upper.median = upper,
              Med = mean(c(lower, upper))))

} # END FUN



#' Partition a vector x into groups based on its quantiles
#'
#' @param x the vector to be partitioned
#' @param cutoffs the quantile cutoffs for the partition. Default are the quartiles.
#' @param quantile.nam logical. Shall the cutoff values be included?
quantile_group <- function(x,
                           cutoffs = c(0.25, 0.5, 0.75),
                           quantile.nam = TRUE){
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

  if(quantile.nam){
    colnames(group.mat) <- nam
  } else{
    colnames(group.mat) <- group.nam
  }
  return(structure(group.mat, type = "quantile_group"))
} # FUN



#' Performs generic ML for a given learning technique (with only one split of the data)
#'
#' @param Z a matrix or data frame of covariates
#' @param D a binary vector of treatment status of length
#' @param Y a vector of responses of length
#' @param propensity.scores a vector of propensity scores
#' @param learner The machine learner that shall be used
#' @param M.set main set
#' @param A.set auxiliary set
#' @param Z_CLAN A matrix of variables that shall be considered for the CLAN. If `NULL` (default), then `Z_CLAN = Z`, i.e. CLAN is performed for all variables in `Z`.
#' @param X1_BLP Specifies the design matrix \eqn{X_1} in the BLP regression. See the documentation of \code{\link{setup_X1}} for details.
#' @param X1_GATES Same as \code{X1_BLP}, just for the the GATES regression.
#' @param HT logical. If TRUE, a HT transformation is applied in BLP and GATES. Default is FALSE.
#' @param vcov_BLP Specifies the covariance matrix estimator in the BLP regression. See the documentation of \code{\link{setup_vcov}} for details.
#' @param vcov_GATES Same as \code{vcov_BLP}, just for the GATES regression.
#' @param equal_variances_CLAN logical. If TRUE, the the two within-group variances of the most and least affected group in CLAN are assumed to be equal. Default is FALSE.
#' @param quantile_cutoffs Cutoff points of quantiles that shall be used for GATES grouping
#' @param diff_GATES Specifies the generic targets of GATES. See the documentation of \code{\link{setup_diff}}.
#' @param diff_CLAN Same as \code{diff_GATES}, just for the CLAN generic targets.
#' @param significance_level Significance level. Default is 0.05
#' @param min_variation minimum variation of the predictions before random noise with distribution N(0, var(Y)/20) is added. Default is 1e-05.
#'
#' TODO: instructions on how mlr3 input is supposed to work (needs to be a string!)
#' TODO: comments on CLAN: If there are categorical variables, apply one-hot-encoding to Z_CLAN. The interpretation then becomes: Is there a factor that is overproportionally present in the least or most affected group?
#'
#' @export
GenericML_single <- function(Z, D, Y,
                             propensity.scores,
                             learner = 'mlr3::lrn("cv_glmnet", s = "lambda.min")',
                             M.set, A.set,
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
                            propensity.scores = propensity.scores,
                            learner = learner,
                            M.set = M.set, A.set = A.set,
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
           propensity.scores,
           learner = 'mlr3::lrn("cv_glmnet", s = "lambda.min")',
           M.set, A.set,
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

    # get the proxy baseline estimator for the main sample
    proxy.baseline.obj <- proxy_baseline_NoChecks(
      Z = Z, D = D, Y = Y,
      auxiliary.sample = A.set,
      learner = learner,
      min_variation = min_variation)
    proxy.baseline     <- proxy.baseline.obj$baseline.predictions.main.sample

    # get the proxy estimator of the CATE for the main sample
    proxy.cate.obj <- proxy_CATE_NoChecks(
      Z = Z, D = D, Y = Y,
      auxiliary.sample = A.set,
      learner = learner,
      proxy.baseline.estimates = proxy.baseline.obj$baseline.predictions.full.sample,
      min_variation = min_variation)
    proxy.cate <- proxy.cate.obj$CATE.predictions.main.sample


    ### step 1b: estimate BLP parameters by OLS ----
    blp.obj <- BLP_NoChecks(
      D = D[M.set],
      Y = Y[M.set],
      propensity.scores  = propensity.scores[M.set],
      proxy.baseline     = proxy.baseline,
      proxy.cate         = proxy.cate,
      HT                 = HT,
      X1_control       = setup_X1(),
      vcov_control       = vcov_BLP,
      significance_level = significance_level)



    ### step 1c: estimate GATES parameters by OLS ----
    # group the proxy estimators for the CATE in the main sample by quantiles
    group.membership.main.sample <- quantile_group(proxy.cate,
                                                   cutoffs = quantile_cutoffs,
                                                   quantile.nam = TRUE)

    gates.obj <- GATES_NoChecks(
      D = D[M.set],
      Y = Y[M.set],
      propensity.scores   = propensity.scores[M.set],
      proxy.baseline      = proxy.baseline,
      proxy.cate          = proxy.cate,
      group.membership.main.sample = group.membership.main.sample,
      HT                  = HT,
      X1_control        = setup_X1(),
      vcov_control        = vcov_GATES,
      diff                = diff_GATES,
      significance_level  = significance_level)


    ### step 1d: estimate CLAN parameters in the main sample ----
    clan.obj <- CLAN_NoChecks(
      Z_CLAN.main.sample = Z_CLAN[M.set,,drop = FALSE],
      group.membership.main.sample = group.membership.main.sample,
      equal.group.variances   = equal_variances_CLAN,
      diff                    = diff_CLAN,
      significance_level      = significance_level)


    ### step 1e: get parameters over which we maximize to find the "best" ML method ----
    best.obj <- lambda_parameters(BLP.obj = blp.obj,
                                  GATES.obj = gates.obj,
                                  proxy.cate.main.sample = proxy.cate,
                                  group.membership.main.sample = group.membership.main.sample)

    ### organize output in a list ----
    return(list(BLP = blp.obj,
                GATES = gates.obj,
                CLAN = clan.obj,
                best = best.obj,
                CATE.proxy = proxy.cate.obj,
                baseline.proxy = proxy.baseline,
                group.membership_M.set = group.membership.main.sample))

  } # FUN



#' returns the two parameters that are used to find the best ML method
#'
#' @param BLP.obj an object as returned by CLAN()
#' @param GATES.obj an object as returned by get.BLP.parameters()
#' @param proxy.cate.main.sample Proxy CATE estimators for the main sample
#' @param group.membership.main.sample a logical matrix with _M_ rows that indicate
#' the group memberships (such a matrix is returned by the function quantile_group())
#' @return lambda and lambda.bar parameters
#'
#' @export
lambda_parameters <- function(BLP.obj,
                              GATES.obj,
                              proxy.cate.main.sample,
                              group.membership.main.sample){

  return(list(lambda = as.numeric(BLP.obj$blp.coefficients["beta.2"]^2 * stats::var(proxy.cate.main.sample)),
              lambda.bar = as.numeric(colMeans(group.membership.main.sample) %*%  GATES.obj$gates.coefficients^2)))

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
