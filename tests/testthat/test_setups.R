# context("Input checks for the setup functions")

# load package
library("GenericML", quietly = TRUE)


test_that("Illegal inputs to setup_diff()", {

  # stopifnot errors
  expect_error(setup_diff("foo"))
  expect_error(setup_diff(c("most", "least")))
  expect_error(setup_diff(subtracted = "foo"))
  expect_error(setup_diff(subtracted = as.matrix(1)))

})


test_that("Illegal inputs to setup_X1()", {

  # stopifnot errors
  expect_error(setup_diff(c("most", "least")))
  expect_error(setup_diff(subtracted = "foo"))
  expect_error(setup_diff(subtracted = as.matrix(1)))

  expect_error(setup_X1("foo"))
  expect_error(setup_X1(covariates = c(1,2,3)))
  expect_error(setup_X1(covariates = data.frame()),
               "If supplied, 'covariates' must be a numeric matrix. Did you supply a data frame?")
  expect_error(setup_X1(covariates = as.character(matrix(c(1,2,3)))))
  expect_error(setup_X1(fixed_effects = "foo"),
               "If supplied, 'fixed_effects' must be a numeric vector")
  expect_error(setup_X1(fixed_effects = c(1.1)),
               "All elements in the vector 'fixed_effects' must be integer-valued")
  expect_error(setup_X1(fixed_effects = as.matrix(1)),
               "If supplied, 'fixed_effects' must be a numeric vector")

})


test_that("Illegal inputs to setup_vcov()", {

  # stopifnot errors
  expect_error(setup_vcov(function(...) print(1)))
  expect_error(setup_vcov("foo"))
  expect_error(setup_vcov(c("vcovHC", "vcovHAC")))

})
