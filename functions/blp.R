#' Estimates the BLP parameters based on the main sample M. 
#' 
#' @param D a binary vector of treatment status of length _|M|_
#' @param Y a vector of responses of length _|M|_
#' @param propensity.scores a vector of propensity scores of length _|M|_
#' @param proxy.baseline a vector of proxy baseline estimates of length _M_
#' @param proxy.cate a vector of proxy CATE estimates of length _M_
#' @param HT.transformation logical. If TRUE, a HT transformation is applied (BLP2 in the paper). Default is FALSE.
#' @param X1.variables a list controlling the variables that shall be used in the matrix X1. The first element of the list, functions_of_Z, needs to be a subset of c("S", "B", "p"), where "p" corresponds to the propensity scores (default is "B"). The seconds element, custom_covariates, is an optional matrix/data frame of custom covariates that shall be included in X1 (default is NULL). The third element, fixed_effects, is a vector of integers that indicates group membership of the observations: For each group, a fixed effect will be added (default is NULL for no fixed effects). Note that in the final matrix X1, a constant 1 will be silently included so that the regression model has an intercept.
#' @param vcov.control a list with two elements called 'estimator' and 'arguments'. The argument 'estimator' is a string specifying the covariance matrix estimator to be used; specifies a covariance estimator function in the sandwich package (https://cran.r-project.org/web/packages/sandwich/sandwich.pdf). Recommended estimators are "vcovBS", "vcovCL", "vcovHAC", and "vcovHC". Default is 'vcovHC'. The element 'arguments' is a list of arguments that shall be passed to the function specified in the element 'estimator'. Default leads to the (homoskedastic) ordinary least squares covariance matrix estimator. See the reference manual of the sandwich package for details (https://cran.r-project.org/web/packages/sandwich/vignettes/sandwich.pdf).
#' @param significance.level significance level for construction of confidence intervals
#' @return BLP coefficients with inference statements
#' 
#' @export
BLP <- function(D, Y, 
                propensity.scores, 
                proxy.baseline,
                proxy.cate, 
                HT.transformation  = FALSE,
                X1.variables       = list(functions_of_Z = c("B"),
                                          custom_covariates = NULL,
                                          fixed_effects = NULL),
                vcov.control       = list(estimator = "vcovHC",
                                          arguments = list(type = "const")),
                significance.level = 0.05){
  
  
  
  # input check
  input.checks.X1(X1.variables)
  
  # fit model according to strategy 1 or 2 in the paper
  do.call(what = get(ifelse(HT.transformation, "BLP.HT", "BLP.classic")),
          args = list(D = D, Y = Y, 
                      propensity.scores  = propensity.scores, 
                      proxy.baseline     = proxy.baseline,
                      proxy.cate         = proxy.cate, 
                      X1.variables       = X1.variables,
                      vcov.control       = vcov.control,
                      significance.level = significance.level))

} # FUN


# helper function for case when there is no HT transformation used. Wrapped by function "BLP"
BLP.classic <- function(D, Y, propensity.scores, 
                        proxy.baseline, proxy.cate, 
                        X1.variables = list(functions_of_Z = c("B"),
                                            custom_covariates = NULL,
                                            fixed_effects = NULL),
                        vcov.control       = list(estimator = "vcovHC",
                                                  arguments = list(type = "const")),
                        significance.level = 0.05){
  
  # prepare weights
  weights <- 1 / (propensity.scores * (1 - propensity.scores))
  
  # prepare covariate matrix X
  X <- data.frame(get.df.from.X1.variables(functions.of.Z_mat = cbind(S = proxy.cate, 
                                                                      B = proxy.baseline, 
                                                                      p = propensity.scores),
                                           X1.variables = X1.variables), 
                  beta.1 = D - propensity.scores, 
                  beta.2 = (D - propensity.scores) * (proxy.cate - mean(proxy.cate))) 
  
  # fit weighted linear regression by OLS
  blp.obj <- lm(Y ~., data = data.frame(Y, X), weights = weights)
  
  # get estimate of covariance matrix of the error terms
  vcov <- get.vcov(x              = blp.obj,
                   vcov.control   = vcov.control)

  # extract coefficients
  coefficients <- lmtest::coeftest(blp.obj, vcov. = vcov)
  
  # return
  return(list(lm.obj = blp.obj, 
              blp.coefficients = blp.obj$coefficients[c("beta.1", "beta.2")],
              generic.targets = generic.targets_BLP(coefficients, significance.level = significance.level),
              coefficients = coefficients))
  
} # END FUN



# helper function for case when there is a HT transformation used. Wrapped by function "BLP"
BLP.HT <- function(D, Y, propensity.scores, 
                   proxy.baseline, proxy.cate, 
                   X1.variables = list(functions_of_Z = c("B"),
                                       custom_covariates = NULL,
                                       fixed_effects = NULL),
                   vcov.control       = list(estimator = "vcovHC",
                                             arguments = list(type = "const")),
                   significance.level = 0.05){
  
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
    X1H$fixed.effects <- fixed.effects # append the fixed effects
    
  } # IF
  
  # prepare covariate matrix X
  X <- data.frame(X1H, 
                  beta.2 = proxy.cate - mean(proxy.cate)) 
  
  # fit linear regression by OLS (intercept is beta.1)
  blp.obj <- lm(YH ~., data = data.frame(YH = Y*H, X))
  names(blp.obj$coefficients) <- c("beta.1", names(blp.obj$coefficients)[-1])

  # get estimate of covariance matrix of the error terms
  vcov <- get.vcov(x              = blp.obj,
                   vcov.control   = vcov.control)
  
  # extract coefficients
  coefficients <- lmtest::coeftest(blp.obj, vcov. = vcov)
  
  # return
  return(list(lm.obj = blp.obj, 
              blp.coefficients = blp.obj$coefficients[c("beta.1", "beta.2")],
              generic.targets = generic.targets_BLP(coefficients, significance.level = significance.level),
              coefficients = coefficients))
  
} # END FUN


# helper function to calculate the generic targets of BLP
generic.targets_BLP <- function(coeftest.object, significance.level = 0.05){
  
  # helper that computes generic targets beta.1, beta.2 in BLP
  coefficients.temp <- coeftest.object[c("beta.1", "beta.2"), 1:3]
  colnames(coefficients.temp) <- c("Estimate", "Std. Error", "z value")
  p.right <- pnorm(coefficients.temp[,"z value"], lower.tail = FALSE) # right p-value: Pr(Z>z)
  p.left  <- pnorm(coefficients.temp[,"z value"], lower.tail = TRUE)  # left p-value: Pr(Z<z)
  ci.lo   <- coefficients.temp[,"Estimate"] - qnorm(1-significance.level/2) * coefficients.temp[,"Std. Error"]
  ci.up   <- coefficients.temp[,"Estimate"] + qnorm(1-significance.level/2) * coefficients.temp[,"Std. Error"]
  generic.targets <- cbind(coefficients.temp, ci.lo, ci.up, p.left, p.right)
  colnames(generic.targets) <- c("Estimate", "Std. Error", "z value", 
                                 "CB lower", "CB upper", "Pr(<z)", "Pr(>z)")
  generic.targets <- generic.targets[,c("Estimate", "CB lower", "CB upper", 
                                        "Std. Error", "z value", "Pr(<z)", "Pr(>z)")]
  
  return(generic.targets)
  
} # FUN
