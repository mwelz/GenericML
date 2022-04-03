# load package
library("GenericML", quietly = TRUE)

if (require("ranger")) {

  ## generate data
  set.seed(1)
  n  <- 100                                  # number of observations
  p  <- 3                                    # number of covariates
  D  <- rbinom(n, 1, 0.5)                    # random treatment assignment
  Z  <- matrix(runif(n*p), n, p)             # design matrix
  Y0 <- as.numeric(Z %*% rexp(p) + rnorm(n)) # potential outcome without treatment
  Y1 <- 2 + Y0                               # potential outcome under treatment
  Y  <- ifelse(D == 1, Y1, Y0)               # observed outcome

  ## assume that data is grouped and that we want each group present in each split,
  # so we use stratified sampling
  groups <- data.frame(group1 = rbinom(n, 1, 0.2),
                       group2 = rbinom(n, 1, 0.3))
  stratify <- setup_stratify(indt = groups,
                             group = c("group1", "group2"),
                             size = 0.5)

  ## specify learner, seed, and call GenericML()
  learners   <- "random_forest"
  num_splits <- 2
  seed       <- 123

  ## we don't expect an error here
  expect_error(GenericML(Z = Z, D = D, Y = Y, num_splits = num_splits,
                         learners_GenericML = learners,
                         stratify = stratify,
                         parallel = FALSE), NA)

  ## set the necessary 'indt' argument to missing, should cause an error
  stratify_corruped <- setup_stratify(foo = groups,
                                      group = c("group1", "group2"),
                                      size = 0.5)
  expect_error(GenericML(Z = Z, D = D, Y = Y, num_splits = num_splits,
                         learners_GenericML = learners,
                         stratify = stratify_corruped,
                         parallel = FALSE))

  ## specify unknown group, should cause an error
  stratify_corruped <- setup_stratify(indt = groups,
                                      group = c("foo1", "group1"),
                                      size = 0.5)
  expect_error(GenericML(Z = Z, D = D, Y = Y, num_splits = num_splits,
                         learners_GenericML = learners,
                         stratify = stratify_corruped,
                         parallel = FALSE))

  ## check if the arguments 'keep.rownames' and 'bothSets' are correctly overwritten
  ## within function. This is indicated by absence of error
  stratify <- setup_stratify(indt = groups,
                             group = c("group1", "group2"),
                             size = 0.5,
                             keep.rownames = FALSE, bothSets = TRUE)
  expect_error(GenericML(Z = Z, D = D, Y = Y, num_splits = num_splits,
                         learners_GenericML = learners,
                         stratify = stratify,
                         parallel = FALSE), NA)

}
