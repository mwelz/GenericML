# load package
library("GenericML", quietly = TRUE)

if (require("ranger")) {

  ## generate data
  set.seed(1)
  n  <- 50                                   # number of observations
  p  <- 3                                    # number of covariates
  D  <- rbinom(n, 1, 0.5)                    # random treatment assignment
  Z  <- matrix(runif(n*p), n, p)             # design matrix
  Y0 <- as.numeric(Z %*% rexp(p) + rnorm(n)) # potential outcome without treatment
  Y1 <- 2 + Y0                               # potential outcome under treatment
  Y  <- ifelse(D == 1, Y1, Y0)               # observed outcome
  learners <- "random_forest"                # specify learner

  set.seed(1)
  x_all <- GenericML(Z = Z, D = D, Y = Y, learners_GenericML = learners,
                     num_splits = 10, seed = 1, parallel = FALSE)

  set.seed(1)
  x <- lapply(1:2, function(...) GenericML(Z = Z, D = D, Y = Y,
                                           learners_GenericML = learners,
                                           num_splits = 5, seed = 1, parallel = FALSE))
  x_comb <- GenericML_combine(x)


  test_that("appending GenericML objects is done correctly in combiner function",
  {
    ## check that appending GenericML objects is done correctly
    expect_equal(x_all$generic_targets$random_forest$BLP[,,2],
                 x_comb$generic_targets$random_forest$BLP[,,2])

    expect_equal(x[[1]]$generic_targets$random_forest$BLP[,,2],
                 x_comb$generic_targets$random_forest$BLP[,,2])

    ## do the same for the second object
    # ideally, we would compare with the 6-th dimension in x_all, but these numbers
    # are different due to different random seeding, because the seeding of x_all cannot be traversed to x
    expect_equal(x_comb$generic_targets$random_forest$BLP[,,6],
                 x[[2]]$generic_targets$random_forest$BLP[,,1])

    ## same for GATES, CLAN, and best
    expect_equal(x_all$generic_targets$random_forest$GATES[,,3],
                 x_comb$generic_targets$random_forest$GATES[,,3])
    expect_equal(x_all$generic_targets$random_forest$CLAN$V1[,,3],
                 x_comb$generic_targets$random_forest$CLAN$V1[,,3])
    expect_equal(x_all$generic_targets$random_forest$best[,,3],
                 x_comb$generic_targets$random_forest$best[,,3])
    expect_equal(x[[2]]$generic_targets$random_forest$GATES[,,4],
                 x_comb$generic_targets$random_forest$GATES[,,9])
    expect_equal(x[[2]]$generic_targets$random_forest$CLAN$V1[,,4],
                 x_comb$generic_targets$random_forest$CLAN$V1[,,9])
    expect_equal(x[[2]]$generic_targets$random_forest$best[,,4],
                 x_comb$generic_targets$random_forest$best[,,9])
  }) # TEST

}
