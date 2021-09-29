
# get lower and upper median as in comment 4.2 in Chernozhukov et al. (2021)
Med <- function(x){
  
  # get the empirical CDF of X
  ecdf.x <- ecdf(x)
  
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
quantile.group <- function(x,
                           cutoffs = c(0.25, 0.5, 0.75),
                           quantile.nam = TRUE){
  # cutoffs are the quantile cutoffs (like c(0.25, 0.5, 0.75))
  
  # get quatiles
  q         <- quantile(x, cutoffs)
  q         <- c(-Inf, q, Inf)
  
  # check if breaks are unique: if x exhibits low variation, there might be empty quantile bins, which can cause an error in the cut() function. In this case, we add random noise to x to induce variation. NB: this bug has been spotted and fixed by Lucas Kitzmueller. All credits for this fix go to him!
  if(length(unique(q)) != length(q)){
    
    # specify standard deviation of the noise (x may have zero variation)
    sd <- ifelse(var(x) == 0, 0.001, sqrt(var(x) / 20))
    
    # add noise and updare quantiles
    x <- x + rnorm(length(x), mean = 0, sd = sd)
    q <- quantile(x, cutoffs)
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
  return(group.mat)
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
#' @param X1.variables_BLP a list controlling the variables that shall be used in the matrix X1 for the BLP regression. The first element of the list, functions_of_Z, needs to be a subset of c("S", "B", "p"), where "p" corresponds to the propensity scores (default is "B"). The second element, custom_covariates, is an optional matrix/data frame of custom covariates that shall be included in X1 (default is NULL). The third element, fixed_effects, is a vector of integers that indicates group membership of the observations: For each group, a fixed effect will be added (default is NULL for no fixed effects). Note that in the final matrix X1, a constant 1 will be silently included so that the regression model has an intercept.
#' @param X1.variables_GATES a list controlling the variables that shall be used in the matrix X1 for the GATES regression. The first element of the list, functions_of_Z, needs to be a subset of c("S", "B", "p"), where "p" corresponds to the propensity scores (default is "B"). The second element, custom_covariates, is an optional matrix/data frame of custom covariates that shall be included in X1 (default is NULL). The third element, fixed_effects, is a vector of integers that indicates group membership of the observations: For each group, a fixed effect will be added (default is NULL for no fixed effects). Note that in the final matrix X1, a constant 1 will be silently included if no HT transformation is applied so that the regression model has an intercept.
#' @param HT.transformation logical. If TRUE, a HT transformation is applied in BLP and GATES. Default is FALSE.
#' @param vcov.control_BLP a list with two elements called 'estimator' and 'arguments'. The argument 'estimator' is a string specifying the covariance matrix estimator to be used in the BLP regression; specifies a covariance estimator function in the sandwich package (https://cran.r-project.org/web/packages/sandwich/sandwich.pdf). Recommended estimators are "vcovBS", "vcovCL", "vcovHAC", and "vcovHC". Default is 'vcovHC'. The element 'arguments' is a list of arguments that shall be passed to the function specified in the element 'estimator'. Default leads to the (homoskedastic) ordinary least squares covariance matrix estimator. See the reference manual of the sandwich package for details (https://cran.r-project.org/web/packages/sandwich/vignettes/sandwich.pdf).
#' @param vcov.control_GATES same as vcov.control_BLP, just for GATES regression
#' @param equal.group.variances_CLAN logical. If TRUE, the the two within-group variances of the most and least affected group in CLAN are assumed to be equal. Default is FALSE.
#' @param quantile.cutoffs Cutoff points of quantiles that shall be used for GATES grouping
#' @param differences.control_GATES a list with two elements called 'group.to.subtract.from' and 'groups.to.be.subtracted'. The first element ('group.to.subtract.from') denotes what shall be the base group to subtract from in GATES; either "most" or "least". The second element ('groups.to.be.subtracted') are the groups to be subtracted from 'group.to.subtract.from', which is a subset of {1,...,K}, where K equals the number of groups.
#' @param differences.control_CLAN same as differences.control_GATES, just for CLAN.
#' @param significance.level Significance level. Default is 0.05
#' @param minimum.variation minimum variation of the predictions before random noise with distribution N(0, var(Y)/20) is added. Default is 1e-05.
#' 
#' TODO: instructions on how mlr3 input is supposed to work (needs to be a string!)
#' TODO: comments on CLAN: If there are categorical variables, apply one-hot-encoding to Z_CLAN. The interpretation then becomes: Is there a factor that is overproportionally present in the least or most affected group?
#' 
#' @export
get.generic.ml.for.given.learner <- function(Z, D, Y, 
                                             propensity.scores,
                                             learner = 'mlr3::lrn("cv_glmnet", s = "lambda.min")',
                                             M.set, A.set,
                                             Z_CLAN                     = NULL, 
                                             X1.variables_BLP           = list(functions_of_Z = c("B"),
                                                                               custom_covariates = NULL,
                                                                               fixed_effects = NULL),
                                             X1.variables_GATES         = list(functions_of_Z = c("B"),
                                                                               custom_covariates = NULL,
                                                                               fixed_effects = NULL),
                                             HT.transformation          = FALSE,
                                             vcov.control_BLP           = list(estimator = "vcovHC",
                                                                            arguments = list(type = "const")),
                                             vcov.control_GATES         = list(estimator = "vcovHC",
                                                                            arguments = list(type = "const")),
                                             equal.group.variances_CLAN = FALSE,
                                             quantile.cutoffs           = c(0.25, 0.5, 0.75),
                                             differences.control_GATES  = list(group.to.subtract.from = "most",
                                                                                groups.to.be.subtracted = 1),
                                             differences.control_CLAN   = list(group.to.subtract.from = "most",
                                                                                groups.to.be.subtracted = 1),  
                                             significance.level         = 0.05,
                                             minimum.variation          = 1e-05){
  
  ### step 1: input checks ---- 
  if(is.null(Z_CLAN)) Z_CLAN <- Z # if no input provided, set it equal to Z
  
  ### step 2a: learn proxy predictors by using the auxiliary set ----
  
  # get the proxy baseline estimator for the main sample
  proxy.baseline.obj <- baseline.proxy.estimator(Z = Z, D = D, Y = Y, 
                                                 auxiliary.sample = A.set, 
                                                 learner = make.mlr3.string(learner, regr = TRUE), 
                                                 minimum.variation = minimum.variation)
  proxy.baseline     <- proxy.baseline.obj$baseline.predictions.main.sample
  
  # get the proxy estimator of the CATE for the main sample
  proxy.cate.obj <- 
    CATE.proxy.estimator(Z = Z, D = D, Y = Y,
                         auxiliary.sample = A.set, 
                         learner = make.mlr3.string(learner, regr = TRUE),
                         proxy.baseline.estimates = proxy.baseline.obj$baseline.predictions.full.sample, 
                         minimum.variation = minimum.variation)
  proxy.cate <- proxy.cate.obj$CATE.predictions.main.sample

  
  ### step 2b: estimate BLP parameters by OLS ----
  blp.obj <- BLP(D = D[M.set], 
                 Y = Y[M.set], 
                 propensity.scores  = propensity.scores[M.set], 
                 proxy.baseline     = proxy.baseline,
                 proxy.cate         = proxy.cate, 
                 HT.transformation  = HT.transformation,
                 X1.variables       = list(functions_of_Z = X1.variables_BLP$functions_of_Z,
                                           custom_covariates = X1.variables_BLP$custom_covariates[M.set,],
                                           fixed_effects = X1.variables_BLP$fixed_effects[M.set]),
                 vcov.control       = vcov.control_BLP,
                 significance.level = significance.level)
  

  ### step 2c: estimate GATES parameters by OLS ----
  # group the proxy estimators for the CATE in the main sample by quantiles
  group.membership.main.sample <- quantile.group(proxy.cate, 
                                                 cutoffs = quantile.cutoffs, 
                                                 quantile.nam = TRUE) 
  
  gates.obj <- GATES(D = D[M.set],
                     Y = Y[M.set], 
                     propensity.scores   = propensity.scores[M.set], 
                     proxy.baseline      = proxy.baseline,
                     proxy.cate          = proxy.cate,
                     group.membership.main.sample = group.membership.main.sample,
                     HT.transformation   = HT.transformation,
                     X1.variables        = list(functions_of_Z = X1.variables_GATES$functions_of_Z,
                                               custom_covariates = X1.variables_GATES$custom_covariates[M.set,],
                                               fixed_effects = X1.variables_GATES$fixed_effects[M.set]),
                     vcov.control        = vcov.control_GATES,
                     differences.control = differences.control_GATES,
                     significance.level  = significance.level)
  
  
  ### step 2d: estimate CLAN parameters in the main sample ----
  clan.obj <- CLAN(Z_CLAN.main.sample = Z_CLAN[M.set,,drop = FALSE], 
                   group.membership.main.sample = group.membership.main.sample, 
                   equal.group.variances   = equal.group.variances_CLAN,
                   differences.control     = differences.control_CLAN,
                   significance.level      = significance.level)
  
  
  ### step 2e: get parameters over which we maximize to find the "best" ML method ----
  best.obj <- best.ml.method.parameters(BLP.obj = blp.obj, 
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
  
} # END FUN


#' returns the two parameters that are used to find the best ML method
#' 
#' @param BLP.obj an object as returned by CLAN()
#' @param GATES.obj an object as returned by get.BLP.parameters()
#' @param proxy.cate.main.sample Proxy CATE estimators for the main sample
#' @param group.membership.main.sample a logical matrix with _M_ rows that indicate 
#' the group memberships (such a matrix is returned by the function quantile.group())
#' @return lambda and lambda.bar parameters
#' 
#' @export
best.ml.method.parameters <- function(BLP.obj,
                                      GATES.obj, 
                                      proxy.cate.main.sample, 
                                      group.membership.main.sample){
  
  return(list(lambda = as.numeric(BLP.obj$blp.coefficients["beta.2"]^2 * var(proxy.cate.main.sample)),
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
#' @export
stratified <- function(group, relative.size, select = NULL){
  
  # perform stratified sampling
  smpl <- splitstackshape::stratified(indt = group, group = group, size = relative.size, select = select)
  return(sort(smpl$ID))
  
} # FUN