#' Estimates the BLP parameters based on the main sample M. 
#' 
#' @param D a binary vector of treatment status of length _|M|_
#' @param Y a vector of responses of length _|M|_
#' @param propensity.scores a vector of propensity scores of length _|M|_
#' @param proxy.baseline a vector of proxy baseline estimates of length _M_
#' @param proxy.cate a vector of proxy CATE estimates of length _M_
#' @param HT.transformation logical. If TRUE, a HT transformation is applied (BLP2 in the paper). Default is FALSE.
#' @param X1.variables a character string specifying the variables in the matrix X1. Needs to be a subset of c("S", "B", "p"), where "p" corresponds to the propensity scores. If no HT transformation is applied, a constant 1 is silently included in X1.
#' @param vcov.type a character string specifying the estimation type of the error covariance matrix. See sandwich::vcovHC for details. Default is "const" (for homoskedasticity)
#' @param significance.level significance level for construction of confidence intervals
#' @return BLP coefficients with inference statements
#' 
#' @export
BLP <- function(D, Y, 
                propensity.scores, 
                proxy.baseline,
                proxy.cate, 
                HT.transformation  = FALSE,
                X1.variables       = c("B"),
                vcov.type          = "const",
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
                      vcov.type          = vcov.type,
                      significance.level = significance.level))

} # FUN


# helper function for case when there is no HT transformation used. Wrapped by function "BLP"
BLP.classic <- function(D, Y, propensity.scores, 
                        proxy.baseline, proxy.cate, 
                        X1.variables = c("B"),
                        vcov.type = "const",
                        significance.level = 0.05){
  
  # prepare weights
  weights <- 1 / (propensity.scores * (1 - propensity.scores))
  
  # prepare matrix X1
  X1.big <- cbind(S = proxy.cate, B = proxy.baseline, p = propensity.scores)
  X1     <- X1.big[, X1.variables]

  # prepare covariate matrix X
  X <- data.frame(X1, 
                  beta.1 = D - propensity.scores, 
                  beta.2 = (D - propensity.scores) * (proxy.cate - mean(proxy.cate))) 
  
  # fit weighted linear regression by OLS
  blp.obj <- lm(Y ~., data = data.frame(Y, X), weights = weights)
  
  # get estimate of error covariance matrix (potentially heteroskedasticity robust)
  vcov <- sandwich::vcovHC(x = blp.obj, type = vcov.type)
  
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
                   X1.variables = c("B"),
                   vcov.type = "const",
                   significance.level = 0.05){
  
  # HT transformation
  H <- (D - propensity.scores) / (propensity.scores * (1 - propensity.scores))
  
  # prepare matrix X1
  X1.big <- cbind(S = proxy.cate, B = proxy.baseline, p = propensity.scores)
  X1     <- X1.big[, X1.variables]
  
  # matrix X_1 * H
  X1H <- X1 * H
  colnames(X1H) <- paste0(colnames(X1), "*H")
  
  # prepare covariate matrix X
  X <- data.frame(X1H, 
                  beta.1 = 1,
                  beta.2 = proxy.cate - mean(proxy.cate)) 
  
  # fit linear regression by OLS (intercept is added implicitly by column of ones in X)
  blp.obj <- lm(formula =  as.formula(paste0("YH ~ ", paste0(colnames(X), collapse = " + "), " + 0")),
                data = data.frame(YH = Y*H, X))

  # get estimate of error covariance matrix (potentially heteroskedasticity robust)
  vcov <- sandwich::vcovHC(x = blp.obj, type = vcov.type)
  
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


# helper that throws error in case of illegal input in 'X1.variables
input.checks.X1 <- function(X1.variables){
  
  legalinput <- X1.variables %in% c("S", "B", "p")
  
  if(!all(legalinput)){
    
    stop(paste0("Entries '", 
                X1.variables[!legalinput], 
                "' of 'X1.variables' are not contained in c('S', 'B', 'p')!"))
    
  } # IF
} # FUN
