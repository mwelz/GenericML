if(requireNamespace("testthat", quietly = TRUE)) {
  library("GenericML", quietly = TRUE)
  testthat::test_check("GenericML")
}
