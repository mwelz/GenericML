# GenericML 0.2.0

- Replaced `1:length(x)`-like loops with safer `seq()`-based counterparts.
- Replaced `if()` conditions comparing `class()` to string with the safer `isa()`.
- Parallel computing is now also supported on Windows.
- Added a method `setup_plot()` that returns the data frame that is used for plotting. Also, made the addition of ATEs in plots optional via the argument `ATE` in `plot.GenericML()`.
- Added a function `GenericML_combine`, which combines multiple `GenericML` objects into one.
- Implemented stratified sampling for sample splitting.


# GenericML 0.1.1

- Fixed a few typos in the documentation.
- Added conditions so that learners based on the package `glmnet` in the tests and examples will be skipped on Solaris machines. Note that this does not prevent an error on Solaris because glmnet is still a `Suggest` of `GenericML` and `glmnet` v4.1.3 cannot be reliably installed on Solaris machines.

# GenericML 0.1.0

- Initial release on CRAN (Nov. 24, 2021)
