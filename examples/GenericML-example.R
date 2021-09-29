#' ------------------------------------------------------------
#' Demonstration of our implementation of Generic ML.
#' Note that this implementation is **not** yet an R package, although we intend to create a package for the CRAN based on it.
#' 
#' Author: mwelz & aalfons & mdemirer
#' Last changed: Sep 28, 2021
#' ------------------------------------------------------------
rm(list = ls()) ; gc(); cat("\014")

# Install the required packages if they are not already installed
required.packages <- c("ggplot2", "mlr3", "mlr3learners", 
                       "mvtnorm", "glmnet", "ranger", "e1071",
                       "kknn", "MASS", "nnet", "xgboost", "sandwich", "lmtest")
new.packages      <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
rm(required.packages, new.packages)


# load the functions
source(paste0(getwd(), "/functions/loader.R"))

### 1. Data Generation (linear, no treatment effect heterogeneity) ----
# We generate $n=5,000$ samples that adhere to a simple linear data generating process. We emulate a randomized experiment. There is no treatment effect heterogeneity since the treatment effect is constant at value two. Hence, Generic ML should not indicate the existence of treatment effect heterogeneity.
set.seed(1)
num.obs  <- 5000
num.vars <- 5

# ATE parameter
ATE <- 2

# random treatment assignment 
D <- rbinom(num.obs, 1, 0.5) 

# covariates
Z <- mvtnorm::rmvnorm(num.obs, mean = rep(0, num.vars), sigma = diag(num.vars))
colnames(Z) <- paste0("z", 1:num.vars)

# counterfactual outcomes
Y0 <- as.numeric(Z %*% c(2, -3, 0, 1, 2)) + rnorm(num.obs)
Y1 <- ATE + Y0

# observed outcome
Y  <- ifelse(D == 1, Y1, Y0) 


### 2. Prepare the arguments for genericML() ---- 

# quantile cutoffs for the GATES grouping of the estimated CATEs 
quantile.cutoffs         <- c(0.2, 0.4, 0.6, 0.8) # 20%, 40%, 60%, 80% quantiles

# specify the learner of the propensity score (non-penalized logistic regression here). Propensity scores can also directly be supplied.
learner.propensity.score <- "mlr3::lrn('glmnet', lambda = 0, alpha = 1)" 

# specify the considered learners of the BCA and the CATE (here: elastic net, random forest, and SVM)
learners.genericML       <- c("elastic.net", "mlr3::lrn('ranger', num.trees = 100)", "mlr3::lrn('svm')")

# specify the data that shall be used for the CLAN
# here, we use all variables of Z and uniformly distributed random noise
Z_CLAN <- cbind(Z, random = runif(num.obs))

# specify the number of splits 
num.splits               <- 100

# specify if a HT transformation shall be used when estimating BLP and GATES
HT.transformation <- FALSE

# A list controlling the variables that shall be used in the matrix X1 for the BLP and GATES regressions. The first element of the list, functions_of_Z, needs to be a subset of c("S", "B", "p"), where "p" corresponds to the propensity scores (default is "B"). The seconds element, custom_covariates, is an optional matrix/data frame of custom covariates that shall be included in X1 (default is NULL). The third element, fixed_effects, is a vector of integers, strings, or a factor thereof that indicates group membership of the observations: For each group, a fixed effect will be added (default is NULL). Note that in the final matrix X1, a constant 1 will be silently included so that the regression model has an intercept (unless HT transformation is applied in GATES).
X1.variables_BLP    <- list(functions_of_Z = c("B"),
                            custom_covariates = NULL,
                            fixed_effects = NULL)
X1.variables_GATES  <- list(functions_of_Z = c("B"),
                            custom_covariates = NULL,
                            fixed_effects = NULL)

# consider differences between group K (most affected) with groups 1 and 2, respectively.
differences.control_GATES  <- list(group.to.subtract.from = "most",
                                   groups.to.be.subtracted = c(1,2)) 
differences.control_CLAN  <- list(group.to.subtract.from = "most",
                                  groups.to.be.subtracted = c(1,2))

# specify the significance level
significance.level       <- 0.05

# specify minimum variation of predictions before Gaussian noise with variance var(Y)/20 is added.
minimum.variation <- 1e-05

# specify which estimator of the error covariance matrix shall be used in BLP and GATES (standard OLS covariance matrix estimator here)
vcov.control_BLP   <- list(estimator = "vcovHC",
                        arguments = list(type = "const"))
vcov.control_GATES <- list(estimator = "vcovHC",
                        arguments = list(type = "const"))

# specify whether of not it should be assumed that the group variances of the most and least affected groups are equal in CLAN. 
equal.group.variances_CLAN <- FALSE

# specify the proportion of samples that shall be selected in the main set
proportion.in.main.set   <- 0.5

# specify whether or not the splits and auxiliary results of the learners shall be stored
store.splits             <- FALSE
store.learners           <- FALSE

### 3. run the genericML() functions with these arguments ----
# runtime: ~121 seconds with R version 4.1.0 on a Dell Latitude 5300 (i5-8265U CPU @ 1.60GHz × 8, 32GB RAM), running on Ubuntu 20.10. Returns a GenericML object.
genML <- GenericML(Z = Z, D = D, Y = Y, 
                   learner.propensity.score = learner.propensity.score, 
                   learners.genericML = learners.genericML,
                   num.splits = num.splits,
                   Z_CLAN = Z_CLAN,
                   HT.transformation = HT.transformation, 
                   X1.variables_BLP = X1.variables_BLP, 
                   X1.variables_GATES = X1.variables_GATES,
                   vcov.control_BLP = vcov.control_BLP, 
                   vcov.control_GATES = vcov.control_GATES,
                   quantile.cutoffs = quantile.cutoffs, 
                   differences.control_GATES = differences.control_GATES,
                   differences.control_CLAN = differences.control_CLAN,
                   equal.group.variances_CLAN = equal.group.variances_CLAN,
                   proportion.in.main.set = proportion.in.main.set, 
                   significance.level = significance.level,
                   minimum.variation = minimum.variation,
                   store.splits = store.splits,
                   store.learners = store.learners)

# save relevant objects
save(genML, ATE, Y, D, Z, file = paste0(getwd(), "/examples/GenericML-example.Rdata"))

# if you don't want to run the genericML() function above, uncomment the line below
# load(file = paste0(getwd(), "/examples/GenericML-example.Rdata"))


### 4. analyze the output ----
# the line below returns the medians of the estimated  \Lambda and \bar{\Lambda}
genML$best.learners$lambda.overview
#                                            lambda lambda.bar
# elastic.net                          0.0003603402   3.981820
# mlr3::lrn('ranger', num.trees = 100) 0.0016003558   4.010251
# mlr3::lrn('svm')                     0.0019966309   3.957013


# We can see that the SVM is the best learner for estimating the CATEs, as it maximizes the median of $\hat{\Lambda}$:
genML$best.learners$best.learner.for.CATE
# "mlr3::lrn('svm')"


# Conversely, the random forest is the best learner for the GATES, as it maximizes the median of $\hat{\bar{\Lambda}}$:
genML$best.learners$best.learner.for.GATES
# "mlr3::lrn('ranger', num.trees = 100)"


# VEIN of BLP
round(genML$VEIN$best.learners$BLP, 5)
#        Estimate CB lower CB upper Pr(<z) adjusted Pr(>z) adjusted
# beta.1  1.98654  1.89041  2.08317               1         0.00000
# beta.2  0.01947 -0.14071  0.18628               1         0.78648
# We see that `beta.1` (the estimate of the ATE) is estimated at ~1.99. True ATE is 2, which is contained in the 90% confidence bounds.  Moreover, `beta.2` is clearly not significant (adjusted $p$-values of both one sided tests are much larger than 0.05). Hence, there is (correctly) no indication of treatment effect heterogeneity. Moreover, the function `genericML.plot()` visualizes these results for the BLP:
plot.GenericML(genML, type = "BLP", title = "VEIN of BLP") 


# VEIN of GATES
round(genML$VEIN$best.learners$GATES, 5)
#                 Estimate CB lower CB upper Pr(<z) adjusted Pr(>z) adjusted
# gamma.1          1.99756  1.76517  2.22686               1         0.00000
# gamma.2          1.98240  1.75276  2.21204               1         0.00000
# gamma.3          1.98286  1.74370  2.21676               1         0.00000
# gamma.4          1.99459  1.77094  2.22709               1         0.00000
# gamma.5          2.02006  1.78922  2.24912               1         0.00000
# gamma.K-gamma.1  0.00069 -0.32917  0.32709               1         0.98721
# gamma.5-gamma.2  0.02119 -0.30804  0.34744               1         0.88855
# All point estimates for the $\gamma$ coefficients are close to each other. Difference between most and least affected group is insignificant (adjusted $p$-values of ~1). Hence, there is (correctly) no indication of treatment effect heterogeneity. We again visualize these results with `genericML.plot()`:
plot.GenericML(genML, type = "GATES", title = "VEIN of GATES") 


# VEIN of CLAN for variable 'z1'
genML$VEIN$best.learners$CLAN$z1
#                     Estimate    CB lower   CB upper Pr(<z) adjusted Pr(>z) adjusted
# delta.1          0.007261430 -0.08740091 0.09945609       1.0000000       0.8227289
# delta.2          0.044316241 -0.04544816 0.13238694       1.0000000       0.3263867
# delta.3          0.028375655 -0.05937453 0.11406139       1.0000000       0.5128145
# delta.4          0.004020818 -0.07849039 0.08880387       1.0000000       0.8738944
# delta.5         -0.038853760 -0.12268665 0.04409473       0.3642137       1.0000000
# delta.5-delta.1 -0.032473153 -0.16097983 0.09603352       0.6204056       1.0000000
# delta.5-delta.2 -0.093908691 -0.21812029 0.03004042       0.1383916       1.0000000
# This correctly indicates that there is no heterogeneity along `z1` (all p-values are much larger than 0.05)
plot.GenericML(genML, type = "CLAN", CLAN.variable = "z1", title = "CLAN of 'z1'") 


# VEIN of CLAN for variable 'random'
plot.GenericML(genML, type = "CLAN", CLAN.variable = "random", title = "CLAN of 'random'") 
# Correctly no evidence for heterogeneity along `random`.


### 5. store all plots ----

pdf(file = paste0(getwd(), "/examples/plots/VEIN-BLP.pdf"))
plot.GenericML(genML, type = "BLP", title = "VEIN of BLP") 
dev.off()

pdf(file = paste0(getwd(), "/examples/plots/VEIN-GATES.pdf"))
plot.GenericML(genML, type = "GATES", title = "VEIN of GATES") 
dev.off()

for(varname in colnames(Z_CLAN)){
  pdf(file = paste0(getwd(), "/examples/plots/VEIN-CLAN-", varname, ".pdf"))
  print(plot.GenericML(genML, type = "CLAN", CLAN.variable = varname,
                 title = paste0("CLAN of '", varname, "'"))) 
  dev.off()
}
