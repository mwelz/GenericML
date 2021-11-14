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
  A_set    <- sample(1:n, n/2)               # auxiliary set

  expect_error(proxy_BCA(Z, D, Y, as.matrix(A_set), learners))
  expect_error(proxy_CATE(Z, D, Y, A_set, learners, proxy_BCA = runif(n+1)))
  expect_error(proxy_CATE(Z, D, Y, A_set, learners, proxy_BCA = as.matrix(runif(n))))

}
