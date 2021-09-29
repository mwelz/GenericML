#' Estimates the GATES parameters based on the main sample M. 
#' 
#' @param D a binary vector of treatment status of length _|M|_
#' @param Y a vector of responses of length _|M|_
#' @param propensity.scores a vector of propensity scores of length _|M|_
#' @param proxy.baseline a vector of proxy baseline estimates of length _M_
#' @param proxy.cate a vector of proxy CATE estimates of length _M_
#' @param group.membership.main.sample a logical matrix with _M_ rows that indicate 
#' the group memberships (such a matrix is returned by the function quantile.group())
#' @param HT.transformation logical. If TRUE, a HT transformation is applied (GATES2 in the paper). Default is FALSE.
#' @param X1.variables a list controlling the variables that shall be used in the matrix X1. The first element of the list, functions_of_Z, needs to be a subset of c("S", "B", "p"), where "p" corresponds to the propensity scores (default is "B"). The seconds element, custom_covariates, is an optional matrix/data frame of custom covariates that shall be included in X1 (default is NULL). The third element, fixed_effects, is a vector of integers that indicates group membership of the observations: For each group, a fixed effect will be added (default is NULL for no fixed effects). Note that in the final matrix X1, a constant 1 will be silently included so that the regression model has an intercept.
#' @param vcov.control a list with two elements called 'estimator' and 'arguments'. The argument 'estimator' is a string specifying the covariance matrix estimator to be used; specifies a covariance estimator function in the sandwich package (https://cran.r-project.org/web/packages/sandwich/sandwich.pdf). Recommended estimators are "vcovBS", "vcovCL", "vcovHAC", and "vcovHC". Default is 'vcovHC'. The element 'arguments' is a list of arguments that shall be passed to the function specified in the element 'estimator'. Default leads to the (homoskedastic) ordinary least squares covariance matrix estimator. See the reference manual of the sandwich package for details (https://cran.r-project.org/web/packages/sandwich/vignettes/sandwich.pdf).
#' @param differences.control a list with two elements called 'group.to.subtract.from' and 'groups.to.be.subtracted'. The first element ('group.to.subtract.from') denotes what shall be the base group to subtract from in GATES; either "most" or "least". The second element ('groups.to.be.subtracted') are the groups to be subtracted from 'group.to.subtract.from', which is a subset of {1,...,K}, where K equals the number of groups.
#' @param significance.level significance level for construction of confidence intervals
#' @return GATES coefficients 
#' 
#' @export
GATES <- function(D, Y, 
                  propensity.scores, 
                  proxy.baseline,
                  proxy.cate,
                  group.membership.main.sample,
                  HT.transformation   = FALSE,
                  X1.variables        = list(functions_of_Z = c("B"),
                                            custom_covariates = NULL,
                                            fixed_effects = NULL),
                  vcov.control        = list(estimator = "vcovHC",
                                            arguments = list(type = "const")),
                  differences.control = list(group.to.subtract.from = "most",
                                              groups.to.be.subtracted = 1),
                  significance.level  = 0.05){
  
  # input check
  InputChecks_D(D)
  InputChecks_Y(Y)
  InputChecks_equal.length2(D, Y)
  InputChecks_equal.length3(propensity.scores, proxy.baseline, proxy.cate)
  InputChecks_equal.length2(Y, propensity.scores)
  InputChecks_X1(X1.variables)
  InputChecks_vcov.control(vcov.control)
  InputChecks_differences.control(differences.control, K = ncol(group.membership.main.sample))
  
  # fit model according to strategy 1 or 2 in the paper
  GATES_NoChecks(D = D, Y = Y,
                 propensity.scores   = propensity.scores, 
                 proxy.baseline      = proxy.baseline,
                 proxy.cate          = proxy.cate, 
                 group.membership.main.sample = group.membership.main.sample,
                 X1.variables        = X1.variables,
                 vcov.control        = vcov.control,
                 differences.control = differences.control,
                 significance.level  = significance.level)
  
} # FUN


# helper function that skips the input checks
GATES_NoChecks <- function(D, Y, 
                           propensity.scores, 
                           proxy.baseline,
                           proxy.cate,
                           group.membership.main.sample,
                           HT.transformation   = FALSE,
                           X1.variables        = list(functions_of_Z = c("B"),
                                                      custom_covariates = NULL,
                                                      fixed_effects = NULL),
                           vcov.control        = list(estimator = "vcovHC",
                                                      arguments = list(type = "const")),
                           differences.control = list(group.to.subtract.from = "most",
                                                      groups.to.be.subtracted = 1),
                           significance.level  = 0.05){
  
  # fit model according to strategy 1 or 2 in the paper
  do.call(what = get(ifelse(HT.transformation, "GATES.HT", "GATES.classic")),
          args = list(D = D, Y = Y,
                      propensity.scores   = propensity.scores, 
                      proxy.baseline      = proxy.baseline,
                      proxy.cate          = proxy.cate, 
                      group.membership.main.sample = group.membership.main.sample,
                      X1.variables        = X1.variables,
                      vcov.control        = vcov.control,
                      differences.control = differences.control,
                      significance.level  = significance.level))
  
} # FUN


# helper function for case when there is no HT transformation used. Wrapped by function "GATES"
GATES.classic <- function(D, Y, 
                          propensity.scores, 
                          proxy.baseline, proxy.cate,
                          group.membership.main.sample,
                          X1.variables = list(functions_of_Z = c("B"),
                                              custom_covariates = NULL,
                                              fixed_effects = NULL),
                          vcov.control       = list(estimator = "vcovHC",
                                                    arguments = list(type = "const")),
                          differences.control = list(group.to.subtract.from = "most",
                                                     groups.to.be.subtracted = 1),
                          significance.level = 0.05){
  
  # make the group membership a binary matrix
  groups <- 1 * group.membership.main.sample
  
  # number of groups
  K <- ncol(groups)
  
  # prepare weights
  weights <- 1 / (propensity.scores * (1 - propensity.scores))
  
  # prepare matrix X1
  X1     <- get.df.from.X1.variables(functions.of.Z_mat = cbind(S = proxy.cate, 
                                                                B = proxy.baseline, 
                                                                p = propensity.scores),
                                     X1.variables = X1.variables)
  
  # prepare covariate matrix
  X <- data.frame(X1, 
                  (D - propensity.scores) * groups)
  colnames(X) <- c(colnames(X1), paste0("gamma.", 1:K))
  
  # fit weighted linear regression by OLS
  gates.obj <- lm(Y ~., data = data.frame(Y, X), weights = weights)
  
  # get estimate of covariance matrix of the error terms
  vcov <- get.vcov(x              = gates.obj,
                   vcov.control   = vcov.control)
  
  # extract the relevant coefficients
  coefficients                 <- lmtest::coeftest(gates.obj, vcov. = vcov)
  gates.coefficients           <- coefficients[paste0("gamma.", 1:K), 1]
  gates.coefficients.quantiles <- colnames(groups)
  names(gates.coefficients.quantiles) <- paste0("gamma.", 1:K)
  
  # return
  return(list(lm.obj = gates.obj, 
              gates.coefficients = gates.coefficients,
              gates.coefficients.quantiles = gates.coefficients.quantiles,
              generic.targets = generic.targets_GATES(coeftest.object = coefficients, 
                                                      K = K, 
                                                      vcov = vcov, 
                                                      significance.level = significance.level,
                                                      differences.control = differences.control),
              coefficients = coefficients))
  
} # END FUN


# helper function for case when there is a HT transformation used. Wrapped by function "GATES"
GATES.HT <- function(D, Y, 
                     propensity.scores, 
                     proxy.baseline, proxy.cate,
                     group.membership.main.sample,
                     X1.variables = list(functions_of_Z = c("B"),
                                         custom_covariates = NULL,
                                         fixed_effects = NULL),
                     vcov.control       = list(estimator = "vcovHC",
                                               arguments = list(type = "const")),
                     differences.control = list(group.to.subtract.from = "most",
                                                groups.to.be.subtracted = 1),
                     significance.level = 0.05){
 
  # make the group membership a binary matrix
  groups <- 1 * group.membership.main.sample
  
  # number of groups
  K <- ncol(groups)
  
  # HT transformation
  H <- (D - propensity.scores) / (propensity.scores * (1 - propensity.scores))
  
  # prepare matrix X1
  X1 <- get.df.from.X1.variables(functions.of.Z_mat = cbind(S = proxy.cate, 
                                                            B = proxy.baseline, 
                                                            p = propensity.scores),
                                 X1.variables = X1.variables)
  
  # construct the matrix X1H (the fixed effects are not multiplied by H, if applicable)
  if(is.null(X1.variables$fixed_effects)){
    
    # matrix X_1 * H
    X1H           <- X1 * H
    colnames(X1H) <- paste0(colnames(X1), ".H")
    
  } else{
    
    fixed.effects     <- X1$fixed.effects  # retain the fixed effects 
    X1                <- X1[,!colnames(X1) %in% "fixed.effects", drop = FALSE]
    X1H               <- X1 * H
    colnames(X1H)     <- paste0(colnames(X1), ".H")
    
    # perform one-hot encoding on fixed effect categories and drop first category to avoid dummy trap
    X1H <- data.frame(X1H, 
                      as.matrix(model.matrix(~ fixed.effects + 0, data = fixed.effects))[,-1])
    
  } # IF
  

  # prepare covariate matrix
  X <- data.frame(X1H,  groups)
  colnames(X) <- c(colnames(X1H), paste0("gamma.", 1:K))
  
  # fit linear regression by OLS (no intercept!)
  gates.obj <- lm(formula =  as.formula(paste0("YH ~ ", paste0(colnames(X), collapse = " + "), " + 0")),
                data = data.frame(YH = Y*H, X))
  
  # get estimate of covariance matrix of the error terms
  vcov <- get.vcov(x              = gates.obj,
                   vcov.control   = vcov.control)
  
  # extract the relevant coefficients
  coefficients                 <- lmtest::coeftest(gates.obj, vcov. = vcov)
  gates.coefficients           <- coefficients[paste0("gamma.", 1:K), 1]
  gates.coefficients.quantiles <- colnames(groups)
  names(gates.coefficients.quantiles) <- paste0("gamma.", 1:K)
  
  # return
  return(list(lm.obj = gates.obj, 
              gates.coefficients = gates.coefficients,
              gates.coefficients.quantiles = gates.coefficients.quantiles,
              generic.targets = generic.targets_GATES(coeftest.object = coefficients, 
                                                      K = K, 
                                                      vcov = vcov, 
                                                      significance.level = significance.level,
                                                      differences.control = differences.control),
              coefficients = coefficients))
  
} # END FUN


# helper function to calculate the generic targets of BLP
generic.targets_GATES <- function(coeftest.object, K, vcov,
                                  significance.level = 0.05,
                                  differences.control = differences.control){
  
  # extract controls
  group.to.subtract.from  <- differences.control$group.to.subtract.from
  groups.to.be.subtracted <- differences.control$groups.to.be.subtracted
  
  # extract coefficients
  coefficients.temp <- coeftest.object[paste0("gamma.", 1:K), 1:3]
  colnames(coefficients.temp) <- c("Estimate", "Std. Error", "z value")
  
  # get quantile
  z          <- qnorm(1-significance.level/2)

  # compute relevant statistics
  p.right <- pnorm(coefficients.temp[,"z value"], lower.tail = FALSE) # right p-value: Pr(Z>z)
  p.left  <- pnorm(coefficients.temp[,"z value"], lower.tail = TRUE)  # left p-value: Pr(Z<z)
  ci.lo   <- coefficients.temp[,"Estimate"] - z * coefficients.temp[,"Std. Error"]
  ci.up   <- coefficients.temp[,"Estimate"] + z * coefficients.temp[,"Std. Error"]
  
  # prepare generic targets of the gammas
  generic.targets <- cbind(coefficients.temp, ci.lo, ci.up, p.left, p.right)
  colnames(generic.targets) <- c("Estimate", "Std. Error", "z value", 
                                 "CB lower", "CB upper", "Pr(<z)", "Pr(>z)")
  generic.targets <- generic.targets[,c("Estimate", "CB lower", "CB upper", 
                                        "Std. Error", "z value", "Pr(<z)", "Pr(>z)")]
  
  # GATES differences: estimate mean and variance of a difference
  if(group.to.subtract.from == "least"){
    
    diff <- coeftest.object["gamma.1", "Estimate"] -
      coeftest.object[paste0("gamma.", groups.to.be.subtracted), "Estimate"]
    
    diff.se <- as.numeric(
      sqrt(vcov["gamma.1", "gamma.1"] + 
                      diag(vcov[paste0("gamma.", groups.to.be.subtracted), 
                                paste0("gamma.", groups.to.be.subtracted), drop = FALSE]) -
                      2 * vcov["gamma.1", paste0("gamma.", groups.to.be.subtracted)]))
    
    nam <- paste0("gamma.1-gamma.", groups.to.be.subtracted)
    
  } else{
    
    diff <- coeftest.object[paste0("gamma.", K), "Estimate"] -
      coeftest.object[paste0("gamma.", groups.to.be.subtracted), "Estimate"]
    
    diff.se <- as.numeric(
      sqrt(vcov[paste0("gamma.", K), paste0("gamma.", K)] + 
                      diag(vcov[paste0("gamma.", groups.to.be.subtracted), 
                                paste0("gamma.", groups.to.be.subtracted), drop = FALSE]) -
                      2 * vcov[paste0("gamma.", K), paste0("gamma.", groups.to.be.subtracted)]))
    
    nam <- paste0("gamma.", K, "-gamma.", groups.to.be.subtracted)
    
  } # IF
  
  ci.lo   <- diff - z * diff.se
  ci.up   <- diff + z * diff.se
  zstat   <- diff / diff.se
  p.right <- pnorm(zstat, lower.tail = FALSE) # right p-value: Pr(Z>z)
  p.left  <- pnorm(zstat, lower.tail = TRUE)  # left p-value: Pr(Z<z)
  
  diff.mat <- cbind(diff, ci.lo, ci.up, diff.se, zstat, p.left, p.right)
  rownames(diff.mat) <- nam
  
  # return final matrix
  rbind(generic.targets, diff.mat)
  
} # FUN
