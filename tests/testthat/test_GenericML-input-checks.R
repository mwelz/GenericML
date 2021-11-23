# context("Input checks for GenericML()")

# load package
library("GenericML", quietly = TRUE)

if (require("glmnet") && require("ranger") && require("rpart")) {

## generate data
set.seed(1)
n  <- 50                                   # number of observations
p  <- 3                                    # number of covariates
D  <- rbinom(n, 1, 0.5)                    # random treatment assignment
Z  <- matrix(runif(n*p), n, p)             # design matrix
Y0 <- as.numeric(Z %*% rexp(p) + rnorm(n)) # potential outcome without treatment
Y1 <- 2 + Y0                               # potential outcome under treatment
Y  <- ifelse(D == 1, Y1, Y0)               # observed outcome
Z_CLAN <- Z                                # data for CLAN
learners <- "lasso"                  # specify learner



## 1. illegal data type of Z, D, Y, Z_CLAN ----
test_that("Errors if incorrect data type of Z, D, Y, Z_CLAN", {

  expect_error(GenericML(data.frame(Z), D, Y, learners, parallel = FALSE),
               "Z must be a numeric matrix. Did you supply a data frame?")
  expect_error(GenericML(as.character(Z), D, Y, learners, parallel = FALSE),
               "Z must be a numeric matrix. Did you supply a data frame?")
  expect_error(GenericML(Z, D, Y, Z_CLAN = data.frame(Z_CLAN), learners, parallel = FALSE),
               "Z_CLAN must be a numeric matrix or NULL. Did you supply a data frame?")
  expect_error(GenericML(Z, as.matrix(D), Y, learners, parallel = FALSE),
               "D must be a numeric vector")
  D_nonbinary <- D; D_nonbinary[1] <- 2
  expect_error(GenericML(Z, D_nonbinary, Y, learners, parallel = FALSE),
               "Treatment assignment D is not binary")
  expect_error(GenericML(Z, D, as.matrix(Y), learners, parallel = FALSE),
               "Y must be a numeric vector")

  # incorrect dimension
  expect_error(GenericML(Z[-1,], D, Y, learners, parallel = FALSE))
  expect_error(GenericML(Z, D, Y, learners, Z_CLAN = Z_CLAN[-1,], parallel = FALSE))
  expect_error(GenericML(Z, D[-1], Y, learners, parallel = FALSE))
  expect_error(GenericML(Z, D, Y[-1], learners, parallel = FALSE))


  # simulate missing values
  Z_NA <- Z
  Z_NA[1,1] <- NA_real_
  D_NA <- D
  D_NA[1] <- NA_real_
  Y_NA <- Y
  Y_NA[1] <- NA_real_
  Z_CLAN_NA <- Z_NA
  expect_error(GenericML(Z_NA, D, Y, learners, parallel = FALSE),
               "Z contains missing values")
  expect_error(GenericML(Z, D_NA, Y, learners, parallel = FALSE),
               "D contains missing values")
  expect_error(GenericML(Z, D, Y_NA, learners, parallel = FALSE),
               "Y contains missing values")
  expect_error(GenericML(Z, D, Y, Z_CLAN = Z_CLAN_NA, learners, parallel = FALSE),
               "Z_CLAN contains missing values")

}) # TEST


## 2. illegal input for learners ----
test_that("Errors in machine learner specification", {

  ## if error message not specified, it is a stopifnot() error
  # GenericML learner
  expect_error(GenericML(Z, D, Y, learners_GenericML = c(1,2,3), parallel = FALSE))
  expect_error(GenericML(Z, D, Y, learners_GenericML = environment(), parallel = FALSE))
  expect_error(GenericML(Z, D, Y, learners, parallel = FALSE,
                         learners_GenericML = "mlr3::foo"),
               "'fooregr.' is not an exported object from 'namespace:mlr3'") # illegal mlr3 syntax
  # illegal mlr3 syntax (type of learning procedure specified)
  expect_error(GenericML(Z, D, Y, parallel = FALSE,
                         learners_GenericML = "mlr3::lrn('regr.ranger', num.trees = 100)"),
               "Element with key 'regr.regr.ranger' not found in DictionaryLearner!")


  # propensity score learner
  expect_error(GenericML(Z, D, Y, learners, learner_propensity_score = environment(), parallel = FALSE))
  expect_error(GenericML(Z, D, Y, learners, parallel = FALSE,
                         learner_propensity_score = c("lasso", "random_forest"))) # multiple learners not allowed
  expect_error(GenericML(Z, D, Y, learners, parallel = FALSE,
                         learner_propensity_score = c(1,2,3)),
               "propensity_scores, Y need to have an equal number of observations")
  expect_error(GenericML(Z, D, Y, learners, parallel = FALSE,
                         learner_propensity_score = "foo"))
  expect_error(GenericML(Z, D, Y, learners, parallel = FALSE,
                         learner_propensity_score = "mlr3::foo"),
               "'fooclassif.' is not an exported object from 'namespace:mlr3'") # illegal mlr3 syntax
  # illegal mlr3 syntax (type of learning procedure specified)
  expect_error(GenericML(Z, D, Y, learners, parallel = FALSE,
                         learner_propensity_score = "mlr3::lrn('classif.ranger', num.trees = 100)"),
               "Element with key 'classif.classif.ranger' not found in DictionaryLearner!")

  # no error with correct specification (skip for computational reasons)
  # expect_error(GenericML(Z, D, Y, num_splits = 2, parallel = FALSE,
  #                       learners_GenericML = c("tree", "lasso", "random_forest")), NA)

}) # TEST


## 3. illegal input for the remaining arguments ----
test_that("Errors in remaining arguments", {

  # single splits are not allowed
  expect_error(GenericML(Z, D, Y, learners, num_splits = 1, parallel = FALSE))

  # cutoffs must be in (0, 1)
  expect_error(GenericML(Z, D, Y, learners, parallel = FALSE,
                         quantile_cutoffs = c(0, 2)))

  expect_error(GenericML(Z, D, Y, learners, X1_BLP = list(), parallel = FALSE),
               "X1_BLP must be an instance of setup_X1()")

  expect_error(GenericML(Z, D, Y, learners, X1_GATES = list(), parallel = FALSE),
               "X1_GATES must be an instance of setup_X1()")

  expect_error(GenericML(Z, D, Y, learners, diff_GATES = list(), parallel = FALSE),
               "diff_GATES must be an instance of setup_diff()")

  expect_error(GenericML(Z, D, Y, learners, diff_CLAN = list(), parallel = FALSE),
               "diff_CLAN must be an instance of setup_diff()")

  expect_error(GenericML(Z, D, Y, learners, vcov_BLP = list(), parallel = FALSE),
               "vcov_BLP must be an instance of setup_vcov()")

  expect_error(GenericML(Z, D, Y, learners, vcov_GATES = list(), parallel = FALSE),
               "vcov_GATES must be an instance of setup_vcov()")

  # sanity check: input for the setup functions
  expect_error(GenericML(Z, D, Y, "random_forest", num_splits = 2, parallel = FALSE,
                         X1_BLP = setup_X1(funs_Z = c("p", "S", "B"))), NA)


  # illegal input to diff arguments
  expect_error(GenericML(Z, D, Y, learners, num_splits = 2, parallel = FALSE,
                         diff_GATES = setup_diff("most", 4)),
               "The most affected group cannot be subtracted from itself")

  expect_error(GenericML(Z, D, Y, learners, num_splits = 2, parallel = FALSE,
                         quantile_cutoffs = c(0.1, 0.2, 0.3, 0.6),
                         diff_GATES = setup_diff("most", 5)),
               "The most affected group cannot be subtracted from itself")

  expect_error(GenericML(Z, D, Y, learners, num_splits = 2, parallel = FALSE,
                         diff_GATES = setup_diff("most", 5)))

  expect_error(GenericML(Z, D, Y, learners, num_splits = 2, parallel = FALSE,
                         diff_GATES = setup_diff("most", 0)))

  expect_error(GenericML(Z, D, Y, learners, num_splits = 2, parallel = FALSE,
                         quantile_cutoffs = c(0.1, 0.2, 0.3, 0.6),
                         diff_GATES = setup_diff("most", 6)))

  expect_error(GenericML(Z, D, Y, learners, num_splits = 2, parallel = FALSE,
                         quantile_cutoffs = c(0.1, 0.2, 0.3, 0.6),
                         diff_GATES = setup_diff("least", 6)))

  expect_error(GenericML(Z, D, Y, learners, num_splits = 2, parallel = FALSE,
                         quantile_cutoffs = c(0.1, 0.2, 0.3, 0.6),
                         diff_GATES = setup_diff("least", 1)),
               "The least affected group cannot be subtracted from itself")

  # illegal input to setup_X1
  expect_error(GenericML(Z, D, Y, learners, num_splits = 2, parallel = FALSE,
                         X1_BLP = setup_X1(covariates = Z[-1,])))
  expect_error(GenericML(Z, D, Y, learners, num_splits = 2, parallel = FALSE,
                         X1_BLP = setup_X1(fixed_effects = rep(1, n-1))))

})


## 4. low-signal input
test_that("function should be able to deal with low-signal input",{

  # illegal propensity scores
  expect_error(GenericML(Z, D, Y, "random_forest", parallel = FALSE,
                         learner_propensity_score = rep(1, n)))

  # illegal treatment assignment (throws both warning and error, as it should)
  #D_ill <- rep(1,n)
  #D_ill[1] <- 0
  #expect_error(GenericML(Z, D_ill, Y, "random_forest"))

  # a constant variable
  const <- rep(1, n)
  expect_error(GenericML(Z = cbind(Z, const), D, Y, "random_forest",
                         num_splits = 2, parallel = FALSE), NA)

  # binary data
  Z <- matrix(sample(c(0,1), n*p, replace = TRUE), n, p)
  Y <- sample(c(0,1), n, replace = TRUE)
  expect_error(GenericML(Z = Z, D, Y, "random_forest",
                         num_splits = 2, parallel = FALSE), NA)



})


} # IF
