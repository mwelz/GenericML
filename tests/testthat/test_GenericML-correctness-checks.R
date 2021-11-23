# skip("commputing-intensive sanity check. Takes too long for R CMD check")
# skip_on_cran("commputing-intensive sanity check. Takes too long for R CMD check")
# # this script contains computing-intensive sanity checks that shall not be run by R CMD check
# # the intensity stems from a large number of splits (1000) to rule out seed-dependence
#
# # load package
# library("GenericML", quietly = TRUE)
#
# if (require("glmnet") && require("ranger") && require("e1071")) {
#
#   test_that("Sanity check: robustness against false positives", {
#
#     ## generate data
#     # observe that there is no treatment effect heterogeneity here
#     set.seed(20211114)
#     n  <- 1000                      # number of observations
#     p  <- 3                         # number of covariates
#     D  <- rbinom(n, 1, 0.5)         # random treatment assignment
#     Z  <- matrix(runif(n*p), n, p)  # design matrix
#     colnames(Z) <- paste0("z", 1:p) # column names
#     ATE <- 2                        # average treatment effect
#     Y0 <- runif(n)                  # potential outcome without treatment
#     Y1 <- ATE + Y0                  # potential outcome under treatment
#     Y  <- ifelse(D == 1, Y1, Y0)    # observed outcome
#     sig_level <- 0.05               # significance level
#
#     # specify learners
#     learners <- c("lasso",
#                   "mlr3::lrn('ranger', num.trees = 100)",
#                   "mlr3::lrn('svm')")
#
#     # perform GenericML for many splits to rule out seed-dependence
#     x <- GenericML(Z, D, Y, learners, num_splits = 1000, significance_level = sig_level,
#                    diff_GATES = setup_diff(subtracted = 1:3),
#                    diff_CLAN = setup_diff(subtracted = 1:3))
#
#     # we expect coverage of ATE
#     beta1_cb <- get_BLP(x, plot = FALSE)["beta.1", c(2,3)]
#     expect_true(beta1_cb[1] <= ATE & ATE <= beta1_cb[2])
#
#     # there is no heterogeneity, so beta.2 shouldn't be significant
#     expect_true(!(get_BLP(x, plot = FALSE)["beta.2", 4] <= sig_level))
#
#     # there is no heterogeneity, so we expect no GATES group difference to be significant
#     gates_diff_pvals <- get_GATES(x, plot = FALSE)[5:7, 4]
#     expect_true(all(gates_diff_pvals > sig_level))
#
#     # there is no heterogeneity, so we expect no CLAN group difference to be significant
#     clan_1_diff_pvals <- get_CLAN(x, plot = FALSE, variable = "z1")[5:7, 4]
#     clan_2_diff_pvals <- get_CLAN(x, plot = FALSE, variable = "z2")[5:7, 4]
#     clan_3_diff_pvals <- get_CLAN(x, plot = FALSE, variable = "z3")[5:7, 4]
#     expect_true(all(clan_1_diff_pvals > sig_level))
#     expect_true(all(clan_2_diff_pvals > sig_level))
#     expect_true(all(clan_3_diff_pvals > sig_level))
#
#   }) # TEST
#
#
#   test_that("Sanity check: ability to identify heterogeneity", {
#
#     ## generate data
#     # observe that there is treatment effect heterogeneity
#     set.seed(20211114)
#     n  <- 1000                      # number of observations
#     p  <- 3                         # number of covariates
#     D  <- rbinom(n, 1, 0.5)         # random treatment assignment
#     Z  <- matrix(runif(n*p), n, p)  # design matrix
#     colnames(Z) <- paste0("z", 1:p) # column names
#     Y0 <- as.numeric(Z %*% rexp(p)) # potential outcome without treatment
#     HTE <- 2 * Z[,1] + ifelse(Z[,2] >= 0.5, 1, -1) # HTE
#     ATE <- mean(HTE)                # ATE
#     Y1 <- HTE + Y0                  # potential outcome under treatment
#     Y  <- ifelse(D == 1, Y1, Y0)    # observed outcome
#     sig_level <- 0.05               # significance level
#
#
#     # specify learners
#     learners <- c("lasso",
#                   "mlr3::lrn('ranger', num.trees = 100)",
#                   "mlr3::lrn('svm')")
#
#     # perform GenericML for many splits to rule out seed-dependence
#     x <- GenericML(Z, D, Y, learners, num_splits = 1000, significance_level = sig_level,
#                    diff_GATES = setup_diff(subtracted = 1:3),
#                    diff_CLAN = setup_diff(subtracted = 1:3))
#
#     # we expect coverage of ATE
#     beta1_cb <- get_BLP(x, plot = FALSE)["beta.1", c(2,3)]
#     expect_true(beta1_cb[1] <= ATE & ATE <= beta1_cb[2])
#
#     # there is heterogeneity, so beta.2 should be significant
#     expect_true(get_BLP(x, plot = FALSE)["beta.2", 4] < sig_level)
#
#     # there is heterogeneity, so we expect all GATES group differences to be significant
#     gates_diff_pvals <- get_GATES(x, plot = FALSE)[5:7, 4]
#     expect_true(all(gates_diff_pvals < sig_level))
#
#
#     ## for evaluating the heterogeneity of CLAN, we need a visual inspection
#
#     # there are two 'streams' of heterogeneity along z1, so we expect a diagonal pattern of the groups,
#     # that is, G1 < G3 < G2 < G4
#     plot(y = Z[,1], x = HTE)
#     get_CLAN(x, plot = TRUE, variable = "z1") # indeed
#
#     # the 2 most affected groups should have a much stronger value of z2 than the 2 least affected groups
#     plot(y = Z[,2], x = HTE)
#     get_CLAN(x, plot = TRUE, variable = "z2") # indeed
#
#     # there is no heterogeneity pattern along z3, so all groups should have roughly equal value
#     plot(y = Z[,3], x = HTE)
#     get_CLAN(x, plot = TRUE, variable = "z3") # indeed
#
#
#   })
#
# }
