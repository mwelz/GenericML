# GenericML
R implementation of [Generic Machine Learning (Chernozhukov, V., Demirer, M., Duflo, E., &amp; Fernández-Val, I., 2021)](https://arxiv.org/abs/1712.04802) using the `mlr3` framework. We intend to extend this implementation to a fully-fledged R package for the CRAN. Please note that this implementation is still work in progress and has not yet been thoroughy tested, so we cannot yet guarantee correctness. If you find a bug, please open an issue or let us know via email.

## Installation
```
# install.packages("devtools")
devtools::install_github("mwelz/GenericML")
```

## Example
We generate `n=5,000` samples that adhere to a simple linear data generating process. We emulate a randomized experiment. There is no treatment effect heterogeneity since the treatment effect is constant at value two. Hence, Generic ML should not indicate the existence of treatment effect heterogeneity.

```{r example}

### 1. Data Generation (linear, no treatment effect heterogeneity) ----
library(GenericML)

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


### 2. Prepare the arguments for GenericML() ----

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

# A list controlling the variables that shall be used in the matrix X1 for the BLP and GATES regressions. 
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

# parallelization options (currently only supported on Unix systems)
parallel  <- TRUE
num.cores <- parallel::detectCores() # maximum number
seed      <- 12345
# Note that the number of cores influences the random number stream. Thus, different choices of `num.cores` may lead to different results.



### 3. Run the GenericML() functions with these arguments ----
# runtime: ~30 seconds with R version 4.1.0 on a Dell Latitude 5300 (i5-8265U CPU @ 1.60GHz × 8, 32GB RAM), running on Ubuntu 21.04. Returns a GenericML object.
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
                   parallel = parallel,
                   num.cores = num.cores,
                   seed = seed,
                   store.splits = store.splits,
                   store.learners = store.learners)

### 4. Analyze the output ----
# the line below returns the medians of the estimated  \Lambda and \bar{\Lambda}
genML$best.learners$lambda.overview

# Get best learner for CATE
genML$best.learners$best.learner.for.CATE

# Get best learner for GATES
genML$best.learners$best.learner.for.GATES

# VEIN of BLP
genML$VEIN$best.learners$BLP
plot(genML, type = "BLP", title = "VEIN of BLP") # plot.GenericML() method
# No indication of treatment effect heterogeneity; see beta.2 coefficient

# VEIN of GATES
genML$VEIN$best.learners$GATES
plot(genML, type = "GATES", title = "VEIN of GATES")
# No indication of heterogeneity

# VEIN of CLAN for variable 'z1'
genML$VEIN$best.learners$CLAN$z1
plot(genML, type = "CLAN", CLAN.variable = "z1", title = "CLAN of 'z1'")
# No indication of heterogeneity

```

## TODO
- [ ] Add optional monotonization of the confindence bounds;
- [ ] Make stratified sampling an argument;
- [ ] Write accessor functions;
- [ ] Make user interface homogeneous and double-check documentation for consistency;
- [ ] Release beta on CRAN.

## Authors
Max Welz (m.welz@erasmusmc.nl), Andreas Alfons (aalfons@ese.eur.nl), Mert Demirer (mdemirer@mit.edu), and Victor Chernozhukov (vchern@mit.edu).
