rm(list = ls()) ; gc(); cat("\014")

# Install the required packages if they are not already installed
required.packages <- c("ggplot2", "mlr3", "mlr3learners", 
                       "mvtnorm", "glmnet", "ranger", "e1071",
                       "kknn", "MASS", "nnet", "xgboost")
new.packages      <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
rm(required.packages, new.packages)


# load the functions
source(paste0(getwd(), "/functions/generic-ml-estimation-funs.R"))
source(paste0(getwd(), "/functions/generic-ml-auxiliary-funs.R"))

### 1. Data Generation (linear, no treatment effect heterogeneity) ----
set.seed(1)
num.obs  <- 1000
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

# quantile cutoffs for the GATES grouping of the estiimated CATEs 
quantile.cutoffs         <- c(0.2, 0.4, 0.6, 0.8) # 20%, 40%, 60%, 80% quantiles

# specify the learner of the propensity score (observe the mlr3 syntax)
# any mlr3 learner can be specified (list: https://mlr3learners.mlr-org.com/)
learner.propensity.score <- "mlr3::lrn('glmnet', lambda = 0, alpha = 1)" # non-penalized logistic regression

# specify the considered learners of the BCA and the CATE.  Any mlr3 learner can be specified (list: https://mlr3learners.mlr-org.com/)
# three learners here: elastic net, random forest, and SVM
learners.genericML       <- c("elastic.net", "mlr3::lrn('ranger', num.trees = 100)", "mlr3::lrn('svm')")

# specify the data that shall be used for the CLAN
# here, we use all variables of Z and uniformly distributed noise
Z.clan <- cbind(Z, random = runif(num.obs))

# specify the number of splits 
num.splits               <- 100

# specify the significance level
significance.level       <- 0.05

# specify the proportion of samples that shall be selected in the main set
proportion.in.main.set   <- 0.5

# specify whether or not the splits and auxiliary results of the learners shall be stored
store.splits             <- FALSE
store.learners           <- FALSE

### 3. run the genericML() functions with these arguments ----
# runtime: ~55 seconds with R version 4.1.0 on a Dell Latitude 5300 (i5-8265U CPU @ 1.60GHz × 8, 32GB RAM), running on Ubuntu 20.10

genML <- genericML(Z = Z, D = D, Y = Y, 
                   learner.propensity.score = learner.propensity.score, 
                   learners.genericML = learners.genericML,
                   num.splits = num.splits,
                   Z.clan = Z.clan,
                   quantile.cutoffs = quantile.cutoffs,
                   proportion.in.main.set = proportion.in.main.set, 
                   significance.level = significance.level,
                   store.splits = store.splits,
                   store.learners = store.learners)

# save relevant objects
save(genML, ATE, Y, D, Z, file = paste0(getwd(), "/examples/GenericML-example.Rdata"))

# if you don't want to run the genericML() function above, uncomment the line below
# load(file = paste0(getwd(), "/examples/GenericML-example.Rdata"))


### 4. analyze the output ----
# the genericML object contains two main lists: 
# 1. "best.learners" for information on finding the best learner,
# 2. "VEIN" for information on the VEIN analysis

## 4.1. "best.learners": performance of the different learners ----
# the line below returns the medians of the estimated  \Lambda and \bar{\Lambda}
genML$best.learners$lambda.overview
# we can see that SVM maximizes both lambda criteria, hence it is the best learner for both CATE and GATES:
genML$best.learners$best.learner.for.CATE
genML$best.learners$best.learner.for.GATES

## 4.2. "VEIN" for the VEIN analysis ----
# we can specify which learner's VEIN we want to return. Choose the best here.
# we can return BLP, GATES, and CLAN objects. All have same structure and contain information on the point estimates, confidence bounds, as well as the raw and adjusted p-values. Since the significance level was set to 5%, the confidence bounds are at 90% confidence level (due to uncertainty from sample splitting).

## 4.2.1. VEIN of BLP ----
genML$VEIN$best.learners$BLP
# beta.1 (the ATE) is estimated at ~1.895. True ATE is 2, which is contained in the 90% CBs.
# beta.2 is clearly not significant (adjusted p-value of ~0.98). Hence, there is (correctly) no indication of treatment effect heterogeneity

## 4.2.2. VEIN of GATES ----
genML$VEIN$best.learners$GATES
# all point estimates for the gamma coefficients are close to each other. Difference between most and least affected group is insignificant (adjusted p-value of ~0.92). Hence, there is (correctly) no indication of treatment effect heterogeneity

## 4.2.3. VEIN of CLAN ----
# VEIN is performed for all variables in the object Z.clan
genML$VEIN$best.learners$CLAN$z1
genML$VEIN$best.learners$CLAN$z2
genML$VEIN$best.learners$CLAN$z3
genML$VEIN$best.learners$CLAN$z4
genML$VEIN$best.learners$CLAN$z5
genML$VEIN$best.learners$CLAN$random
# there does not seem to be heterogeneity along any variable in Z.clan (and rightfully so)

## 5. Visualization of the output ----
library(ggplot2)
# the function genericML.plot() visualizes each of the VEIN analyses
genericML.plot(genML, type = "GATES") # no hetero


# analyze
genML$VEIN$best.learners$GATES # difference is insignificant, so no hetero
genML$VEIN$best.learners$BLP  # beta2 is insignificant, so no hetero
genML$VEIN$best.learners$CLAN$z1 # there seems to be hetero along z1

# GATES
genericML.plot(genML, type = "GATES") # no hetero
genericML.plot(genML, type = "BLP")   # no hetero

genericML.plot(genML, type = "CLAN", CLAN.variable = "z1")   # no hetero
genericML.plot(genML, type = "CLAN", CLAN.variable = "z2")   # hetero
genericML.plot(genML, type = "CLAN", CLAN.variable = "z3")   # slight hetero
genericML.plot(genML, type = "CLAN", CLAN.variable = "z4")   # no hetero
genericML.plot(genML, type = "CLAN", CLAN.variable = "z5")   # no hetero
