# GenericML: Generic Machine Learning Inference

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/GenericML)](https://cran.r-project.org/package=GenericML)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/GenericML)](https://cran.r-project.org/package=GenericML)

To cite `GenericML` in publications, please use:

> Welz M., Alfons, A., Demirer, M. and Chernozhukov, V. (2022). `GenericML`: Generic Machine Learning Inference. `R` package version 0.2.3. URL: https://CRAN.R-project.org/package=GenericML.


## Summary

`R` implementation of Generic Machine Learning Inference on heterogeneous treatment effects in randomized experiments as proposed in [Chernozhukov, Demirer, Duflo and Fernández-Val (2020)](https://arxiv.org/abs/1712.04802). This package's workhorse is the `mlr3` framework of [Lang et al. (2019)](https://joss.theoj.org/papers/10.21105/joss.01903), which enables the specification of a wide variety of machine learners. The main functionality, `GenericML()`, runs Algorithm 1 in [Chernozhukov, Demirer, Duflo and Fernández-Val (2020)](https://arxiv.org/abs/1712.04802) for a suite of user-specified machine learners. All steps in the algorithm are customizable via setup functions. Methods for printing and plotting are available for objects returned by `GenericML()`. Parallel computing is supported.

## Installation

### From CRAN
The package `GenericML` is on the CRAN (The Comprehensive R Archive Network), hence the latest release can be easily installed from the `R` command line via
```R
install.packages("GenericML")
```

### Building from source

To install the latest (possibly unstable) development version from GitHub, you can pull this repository and install it from the `R` command line via
```R
install.packages("devtools")
devtools::install_github("mwelz/GenericML")
```
If you already have the package `devtools` installed, you can skip the first line.

## Community guidelines

### Report issues and request features

If you experience any bugs or issues or if you have any suggestions for additional features, please submit an issue via the [*Issues*](https://github.com/mwelz/GenericML/issues) tab of this repository. Please have a look at existing issues first to see if your problem or feature request has already been discussed.

### Contribute to the package

If you want to contribute to the package, you can fork this repository and create a pull request after implementing the desired functionality.

### Ask for help

If you need help using the package, or if you are interested in collaborations related to this project, please get in touch with the [package maintainer](https://mwelz.github.io/).


## Example
We generate `n = 1000` samples from a randomized experiment which exhibit substantial treatment effect heterogeneity. Hence, Generic ML should indicate the existence of treatment effect heterogeneity and also capture the heterogeneity patterns.

### Data generation and main functionality

```R
library(GenericML)

### 1. Data Generation ----
set.seed(31684591)

n  <- 1000                      # number of observations
p  <- 3                         # number of covariates
D  <- rbinom(n, 1, 0.5)         # random treatment assignment
Z  <- matrix(runif(n*p), n, p)  # design matrix
colnames(Z) <- paste0("z", 1:p) # column names
Y0 <- as.numeric(Z %*% rexp(p)) # potential outcome without treatment

## simulate heterogeneous treatment effect
# treatment effect increases with Z1, has a level shift along Z2 (at 0.5), and has no pattern along Z3
HTE <- 2 * Z[,1] + ifelse(Z[,2] >= 0.5, 1, -1)
ATE <- mean(HTE)                # average treatment effect
Y1  <- HTE + Y0                 # potential outcome under treatment
Y   <- ifelse(D == 1, Y1, Y0)   # observed outcome


### 2. Prepare the arguments for GenericML() ----

# quantile cutoffs for the GATES grouping of the estimated CATEs
quantile_cutoffs <- c(0.25, 0.5, 0.75) # 25%, 50%, and 75% quantiles

# specify the learner of the propensity score (non-penalized logistic regression here). Propensity scores can also directly be supplied.
learner_propensity_score <- "mlr3::lrn('glmnet', lambda = 0, alpha = 1)"

# specify the considered learners of the BCA and the CATE (here: lasso, random forest, and SVM)
learners_GenericML <- c("lasso", "mlr3::lrn('ranger', num.trees = 100)", "mlr3::lrn('svm')")

# specify the data that shall be used for the CLAN
# here, we use all variables of Z and uniformly distributed random noise
Z_CLAN <- cbind(Z, random = runif(n))

# specify the number of splits (many to rule out seed-dependence of results)
num_splits <- 1000

# specify if a HT transformation shall be used when estimating BLP and GATES
HT <- FALSE

# A list controlling the variables that shall be used in the matrix X1 for the BLP and GATES regressions.
X1_BLP   <- setup_X1()
X1_GATES <- setup_X1()

# consider differences between group K (most affected) with groups 1, 2, and 3, respectively.
diff_GATES <- setup_diff(subtract_from = "most",
                         subtracted = 1:3)
diff_CLAN  <- setup_diff(subtract_from = "most",
                         subtracted = 1:3)

# specify the significance level
significance_level <- 0.05

# specify minimum variation of predictions before Gaussian noise with variance var(Y)/20 is added.
min_variation <- 1e-05

# specify which estimator of the error covariance matrix shall be used in BLP and GATES (standard OLS covariance matrix estimator here)
vcov_BLP   <- setup_vcov()
vcov_GATES <- setup_vcov()

# specify whether of not it should be assumed that the group variances of the most and least affected groups are equal in CLAN.
equal_variances_CLAN <- FALSE

# specify the proportion of samples that shall be selected in the auxiliary set
prop_aux <- 0.5

# specify sampling strategy (possibly stratified). Here ordinary random sampling is used.
stratify <- setup_stratify()

# specify whether or not the splits and auxiliary results of the learners shall be stored
store_splits   <- TRUE
store_learners <- FALSE # to save memory

# parallelization options
parallel  <- TRUE
num_cores <- 8      # 8 cores
seed      <- 123456
# Note that the number of cores as well as your type of operating system (Unix vs. Windows) influences the random number stream. Thus, different choices of `num_cores` may lead to different results. Results of parallel processes are reproducible across all Unix systems, but might deviate on Windows systems.



### 3. Run the GenericML() function with these arguments ----
# runtime: ~90 seconds with R version 4.2.0 on a Dell Latitude 5300 (i5-8265U CPU @ 1.60GHz × 8, 32GB RAM), running on Ubuntu 21.10. Returns a GenericML object.
x <- GenericML(Z = Z, D = D, Y = Y,
               learner_propensity_score = learner_propensity_score,
               learners_GenericML = learners_GenericML,
               num_splits = num_splits,
               Z_CLAN = Z_CLAN,
               HT = HT,
               X1_BLP = X1_BLP,
               X1_GATES = X1_GATES,
               vcov_BLP = vcov_BLP,
               vcov_GATES = vcov_GATES,
               quantile_cutoffs = quantile_cutoffs,
               diff_GATES = diff_GATES,
               diff_CLAN = diff_CLAN,
               equal_variances_CLAN = equal_variances_CLAN,
               prop_aux = prop_aux,
               stratify = stratify,
               significance_level = significance_level,
               min_variation = min_variation,
               parallel = parallel,
               num_cores = num_cores,
               seed = seed,
               store_splits = store_splits,
               store_learners = store_learners)
               

### 4. General results ----

## print
genML
# GenericML object with the following specifications:
# 	* Propensity score learner: mlr3::lrn('glmnet', lambda = 0, alpha = 1) 
# 	* Generic ML learners: lasso, mlr3::lrn('ranger', num.trees = 100), mlr3::lrn('svm') 
# 	* S = 1000 splits are used
# 	* No HT transformation is used
# 
# The 90 % confidence intervals of the best BLP estimates are given by
# 	 beta.1: ( 0.894 , 1.056 )	 beta.2: ( 0.944 , 1.095 )
# The best learner for the BLP is mlr3::lrn('svm') (lambda of 1.172)
# The best learner for the GATES and CLAN is mlr3::lrn('ranger', num.trees = 100) (lambda.bar of 2.0818)


## get the medians of the estimated  \Lambda and \bar{\Lambda} to find best learners
get_best(genML)
#                                      lambda lambda.bar
# lasso                                 1.119      2.005
# mlr3::lrn('ranger', num.trees = 100)  1.167      2.082
# mlr3::lrn('svm')                      1.172      2.063
# ---
# The best learner for BLP is mlr3::lrn('svm') with lambda = 1.172.
# The best learner for GATES and CLAN is mlr3::lrn('ranger', num.trees = 100) with lambda.bar = 2.0818.
```

*We emphasize that the number of cores and your type of operating system affect the random number stream. Thus, different choices of `num_cores` may lead to different results. Moreover, results of parallel processes are reproducible across all Unix systems, but might deviate on Windows systems. Consequently, the results below are only reproducible on Unix systems for `num_cores = 8`.*

### Best Linear Predictor (BLP) analysis

We use the `get_BLP()` acceessor function to extract the results of the BLP analysis. We can see from the print and plot that the true ATE of about 0.979 is contained in the 90% confidence interval of `beta.1`. Moreover, we reject the null of no significance of `beta.2` at any reasonable level, which is expected since there is substantial treatment effect heterogeneity.

```R
results_BLP <- get_BLP(genML)
results_BLP # print method
# BLP generic targets
# ---
#        Estimate CI lower CI upper p value
# beta.1   0.9750   0.8939    1.056  <2e-16 ***
# beta.2   1.0198   0.9442    1.095  <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# ---
# Confidence level of confidence interval [CI lower, CI upper]: 90 %

plot(results_BLP) # plot method
```

<img src="./inst/doc/repo_plots/BLP.svg" width="67%" style="display: block; margin: auto;" />

### Group Average Treatment Effects (GATES) analysis

There is treatment effect heterogeneity in the data generating process, so we expect a trend in the GATES per-group estimates as well as significance of all group differences. This is indeed the case.

```R
results_GATES <- get_GATES(genML)
results_GATES # print method
# GATES generic targets
# ---
#                 Estimate CI lower CI upper  p value
# gamma.1          -0.4475  -0.6488   -0.247 1.12e-05 ***
# gamma.2           0.4997   0.2976    0.705 1.18e-06 ***
# gamma.3           1.4308   1.2298    1.634  < 2e-16 ***
# gamma.4           2.4079   2.2063    2.612  < 2e-16 ***
# gamma.4-gamma.1   2.8603   2.5732    3.142  < 2e-16 ***
# gamma.4-gamma.2   1.9056   1.6212    2.192  < 2e-16 ***
# gamma.4-gamma.3   0.9824   0.6923    1.268 1.92e-11 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# ---
# Confidence level of confidence interval [CI lower, CI upper]: 90 %

plot(results_GATES) # plot method
```

<img src="./inst/doc/repo_plots/GATES.svg" width="67%" style="display: block; margin: auto;" />

### Classification Analysis (CLAN) of each covariate

#### CLAN of first covariate

We first inspect the (true) treatment effect heterogeneity along the first variable:

```R
plot(y = Z_CLAN[, "z1"], x = HTE, xlab = "True HTE", ylab = "Value of z1")
```

<img src="./inst/doc/repo_plots/z1.svg" width="67%" style="display: block; margin: auto;" />

The heterogeneity exhibits a jump pattern along the first variable. We thus expect that `G1 < G3`, `G2 < G4`, `G1 = G2`, `G3 = G4`, where the `G` denote the variable’s within-group averages. The groups are formed by treatment effect strength. Let’s see what CLAN suggests:

```R
results_CLAN_z1 <- get_CLAN(genML, variable = "z1")
plot(results_CLAN_z1)
```

<img src="./inst/doc/repo_plots/CLAN_z1.svg" width="67%" style="display: block; margin: auto;" />

CLAN indeed captured the correct pattern and that `G1 < G3`, `G2 < G4`, `G1 = G2`, and `G3 = G4`.

#### CLAN of second covariate

We inspect the (true) treatment effect heterogeneity along the second variable:

```R
plot(y = Z_CLAN[, "z2"], x = HTE, xlab = "True HTE", ylab = "Value of z2")
```

<img src="./inst/doc/repo_plots/z2.svg" width="67%" style="display: block; margin: auto;" />

We clearly see the level shift at (1, 0.5). Thus, we expect that the two most affected groups should have a much stronger value of `z2` than the two least affected groups. Moreover, the two groups `G1` and `G2` should have the same value of `z2` and the two groups `G3` and `G4` should also have the same value. CLAN indeed captures this pattern:

```R
results_CLAN_z2 <- get_CLAN(genML, variable = "z2")
plot(results_CLAN_z2)
```

<img src="./inst/doc/repo_plots/CLAN_z2.svg" width="67%" style="display: block; margin: auto;" />

#### CLAN of third covariate

We inspect the (true) treatment effect heterogeneity along the third variable:

```R
plot(y = Z_CLAN[, "z3"], x = HTE, xlab = "True HTE", ylab = "Value of z3")
```

<img src="./inst/doc/repo_plots/z3.svg" width="67%" style="display: block; margin: auto;" />

There is no heterogeneity pattern along `z3`, so all CLAN groups should have roughly the same value. This is indeed the case:

```R
results_CLAN_z3 <- get_CLAN(genML, variable = "z3")
plot(results_CLAN_z3)
```

<img src="./inst/doc/repo_plots/CLAN_z3.svg" width="67%" style="display: block; margin: auto;" />

#### CLAN of fourth covariate

We inspect the (true) treatment effect heterogeneity along the fourth variable which is just random noise:

```R
plot(y = Z_CLAN[, "random"], x = HTE, xlab = "True HTE", ylab = "Value of 'random'")
```

<img src="./inst/doc/repo_plots/random.svg" width="67%" style="display: block; margin: auto;" />

There is no heterogeneity along `random`, so all CLAN groups should have roughly the same value. This is indeed the case:

```R
results_CLAN_random <- get_CLAN(genML, variable = "random")
plot(results_CLAN_random)
```

<img src="./inst/doc/repo_plots/CLAN_random.svg" width="67%" style="display: block; margin: auto;" />


## Authors
Max Welz (welz@ese.eur.nl), Andreas Alfons (alfons@ese.eur.nl), Mert Demirer (mdemirer@mit.edu), and Victor Chernozhukov (vchern@mit.edu).
