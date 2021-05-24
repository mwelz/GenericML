#' ------------------------------------------------------------
#' Demonstration of our implementation of Generic ML.
#' Note that this implementation is **not** yet an R package, although we intend to create a package for the CRAN based on it.
#' 
#' Author: mwelz & aalfons
#' Last changed: May 23, 2021
#' ------------------------------------------------------------
rm(list = ls()) ; gc(); cat("\014")

# Install the required packages if they are not already installed
required.packages <- c("ggplot2", "mlr3", "mlr3learners", 
                       "mvtnorm", "glmnet", "ranger", "e1071",
                       "kknn", "MASS", "nnet", "xgboost")
new.packages      <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
rm(required.packages, new.packages)
library(ggplot2)


# load the functions
source(paste0(getwd(), "/functions/generic-ml-estimation-funs.R"))
source(paste0(getwd(), "/functions/generic-ml-auxiliary-funs.R"))

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

# specify the learner of the propensity score (non-penalized logistic regression here)
learner.propensity.score <- "mlr3::lrn('glmnet', lambda = 0, alpha = 1)" 

# specify the considered learners of the BCA and the CATE (here: elastic net, random forest, and SVM)
learners.genericML       <- c("elastic.net", "mlr3::lrn('ranger', num.trees = 100)", "mlr3::lrn('svm')")

# specify the data that shall be used for the CLAN
# here, we use all variables of Z and uniformly distributed random noise
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
# runtime: ~121 seconds with R version 4.1.0 on a Dell Latitude 5300 (i5-8265U CPU @ 1.60GHz × 8, 32GB RAM), running on Ubuntu 20.10
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
# the line below returns the medians of the estimated  \Lambda and \bar{\Lambda}
genML$best.learners$lambda.overview
#'                                      lambda         lambda.bar
#' elastic.net                          0.0003662581   9950.862
#' mlr3::lrn('ranger', num.trees = 100) 0.0009451692   9823.487
#' mlr3::lrn('svm')                     0.0011130656   9826.750


# We can see that the SVM is the best learner for estimating the CATEs, as it maximizes the median of $\hat{\Lambda}$:
genML$best.learners$best.learner.for.CATE
#' "mlr3::lrn('svm')"


# Conversely, the elastic net is the best learner for the GATES, as it maximizes the median of $\hat{\bar{\Lambda}}$:
genML$best.learners$best.learner.for.GATES
#' "elastic net"


# VEIN of BLP
round(genML$VEIN$best.learners$BLP, 5)
#'        Estimate CB lower CB upper p-value adjusted p-value raw
#' beta.1  1.97975  1.88779  2.07171          0.00000     0.00000
#' beta.2  0.02613 -0.12937  0.18410          0.94309     0.47154
# We see that `beta.1` (the estimate of the ATE) is estimated at ~1.98. True ATE is 2, which is contained in the 90% confidence bounds.  Moreover, `beta.2` is clearly not significant (adjusted $p$-value of ~0.94). Hence, there is (correctly) no indication of treatment effect heterogeneity. Moreover, the function `genericML.plot()` visualizes these results for the BLP:
genericML.plot(genML, type = "BLP", title = "VEIN of BLP") 


# VEIN of GATES
round(genML$VEIN$best.learners$GATES, 5)
#'                 Estimate CB lower CB upper p-value adjusted p-value raw
#' gamma.1          1.93589  1.75543  2.11599                0     0.00000
#' gamma.2          1.99687  1.81672  2.17643                0     0.00000
#' gamma.3          2.02001  1.84212  2.19790                0     0.00000
#' gamma.4          2.01884  1.83918  2.19775                0     0.00000
#' gamma.5          1.97835  1.79725  2.15980                0     0.00000
#' gamma.K-gamma.1  0.03154 -0.22497  0.28706                1     0.57596
# All point estimates for the $\gamma$ coefficients are close to each other. Difference between most and least affected group is insignificant (adjusted $p$-value of ~1). Hence, there is (correctly) no indication of treatment effect heterogeneity. We again visualize these results with `genericML.plot()`:
genericML.plot(genML, type = "GATES", title = "VEIN of GATES") 


# VEIN of CLAN for variable 'z1'
genML$VEIN$best.learners$CLAN$z1
#'                   Estimate   CB lower   CB upper p-value adjusted  p-value raw
#' delta.1          0.2091063  0.1222465  0.2959660     9.418759e-30 4.709380e-30
#' delta.K         -0.2365256 -0.3227453 -0.1476362     7.966651e-31 3.983326e-31
#' delta.K-delta.1 -0.4673637 -0.5901810 -0.3445465     7.203990e-60 3.601995e-60
# This indicates some evidence for weak heterogeneity along the variable `z1`. This could be due to the fact that `z1` is positively correlated with the outcome Y.
genericML.plot(genML, type = "CLAN", CLAN.variable = "z1", title = "CLAN of 'z1'") 


# VEIN of CLAN for variable 'random'
genericML.plot(genML, type = "CLAN", CLAN.variable = "random", title = "CLAN of 'random'") 
# Correctly no evidence for heterogeneity along `random`.


### 5. store all plots ----

pdf(file = paste0(getwd(), "/examples/plots/VEIN-BLP.pdf"))
genericML.plot(genML, type = "BLP", title = "VEIN of BLP") 
dev.off()

pdf(file = paste0(getwd(), "/examples/plots/VEIN-GATES.pdf"))
genericML.plot(genML, type = "GATES", title = "VEIN of GATES") 
dev.off()

for(varname in colnames(Z.clan)){
  pdf(file = paste0(getwd(), "/examples/plots/VEIN-CLAN-", varname, ".pdf"))
  print(genericML.plot(genML, type = "CLAN", CLAN.variable = varname,
                 title = paste0("CLAN of '", varname, "'"))) 
  dev.off()
}
