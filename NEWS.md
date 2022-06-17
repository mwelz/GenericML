# GenericML 0.2.2

- Added class structure for accessor function objects
- Ensured consistency in documentation.
- Added new function, `heterogeneity_CLAN()`, that investigates the presence of treatment effect heterogeneity along all CLAN variables.
- Added function `get_best()` that returns the best learner.
- Changed behavior of `get_CLAN()` to not plot ATE estimates when `plot = TRUE`.

# GenericML 0.2.1

- Replaced `isa()` with `inherits()` to avoid reliance on `R >= 4.1`.
- Changed default in `parallel` argument in `GenericML` to `FALSE`.

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
