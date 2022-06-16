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

  ## specify learner, number of splits, and seed
  learners   <- "random_forest"
  num_splits <- 2
  seed       <- 123

  ## check 1: reproducibility of serial implementation
  # without stratified sampling
  x1 <- GenericML(Z = Z, D = D, Y = Y, num_splits = num_splits,
                  learners_GenericML = learners,
                  seed = seed,
                  parallel = FALSE)
  x2 <- GenericML(Z = Z, D = D, Y = Y, num_splits = num_splits,
                  learners_GenericML = learners,
                  seed = seed,
                  parallel = FALSE)

  expect_equal(object = x1$VEIN$best_learners$BLP[1,],
               expected = x2$VEIN$best_learners$BLP[1,])

  # with stratified sampling
  x1 <- GenericML(Z = Z, D = D, Y = Y, num_splits = num_splits,
                  learners_GenericML = learners,
                  stratify = stratify, seed = seed,
                  parallel = FALSE)
  x2 <- GenericML(Z = Z, D = D, Y = Y, num_splits = num_splits,
                  learners_GenericML = learners,
                  stratify = stratify, seed = seed,
                  parallel = FALSE)

  expect_equal(object = x1$VEIN$best_learners$BLP[1,],
               expected = x2$VEIN$best_learners$BLP[1,])

  ## on the contrary, no specification of seeds should cause non-reproducible results
  x1 <- GenericML(Z = Z, D = D, Y = Y, num_splits = num_splits,
                  learners_GenericML = learners,
                  stratify = stratify, seed = NULL,
                  parallel = FALSE)
  x2 <- GenericML(Z = Z, D = D, Y = Y, num_splits = num_splits,
                  learners_GenericML = learners,
                  stratify = stratify, seed = NULL,
                  parallel = FALSE)

  expect_false(isTRUE(all.equal(target  = x1$VEIN$best_learners$BLP[1,],
                                current = x2$VEIN$best_learners$BLP[1,])))


  ## check 2: reproducibility of parallel implementation
  ## this is commented out on the CRAN release because CRAN disapproves
  ## tests that require multiple cores
  # num_cores <- 2L
  #
  # # without stratified sampling
  # x1 <- GenericML(Z = Z, D = D, Y = Y, num_splits = num_splits,
  #                 learners_GenericML = learners,
  #                 seed = seed,
  #                 parallel = TRUE, num_cores = num_cores)
  # x2 <- GenericML(Z = Z, D = D, Y = Y, num_splits = num_splits,
  #                 learners_GenericML = learners,
  #                 seed = seed,
  #                 parallel = TRUE, num_cores = num_cores)
  #
  # expect_equal(object = x1$VEIN$best_learners$BLP[1,],
  #              expected = x2$VEIN$best_learners$BLP[1,])
  #
  # # with stratified sampling
  # x1 <- GenericML(Z = Z, D = D, Y = Y, num_splits = num_splits,
  #                 learners_GenericML = learners,
  #                 stratify = stratify, seed = seed,
  #                 parallel = TRUE, num_cores = num_cores)
  # x2 <- GenericML(Z = Z, D = D, Y = Y, num_splits = num_splits,
  #                 learners_GenericML = learners,
  #                 stratify = stratify, seed = seed,
  #                 parallel = TRUE, num_cores = num_cores)
  #
  # expect_equal(object = x1$VEIN$best_learners$BLP[1,],
  #              expected = x2$VEIN$best_learners$BLP[1,])

}
