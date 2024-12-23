library("GenericML")

# load data
load("inst/doc/repo_plots/plot_data.Rdata")

### 4.1 analyze results ----
## print
genML
# GenericML object with the following specifications:
#   * Propensity score learner: mlr3::lrn('glmnet', lambda = 0, alpha = 1)
# * Generic ML learners: lasso, mlr3::lrn('ranger', num.trees = 100), mlr3::lrn('svm')
# * S = 1000 splits are used
# * No HT transformation is used
#
# The 90 % confidence intervals of the best BLP estimates are given by
# beta.1: ( 0.894 , 1.056 )	 beta.2: ( 0.944 , 1.095 )
# The best learner for the BLP is mlr3::lrn('svm') (lambda of 1.172)
# The best learner for the GATES and CLAN is mlr3::lrn('ranger', num.trees = 100) (lambda.bar of 2.0818)


## get the medians of the estimated  \Lambda and \bar{\Lambda} to find best learners
get_best(genML)
# lambda lambda.bar
# lasso                                 1.119      2.005
# mlr3::lrn('ranger', num.trees = 100)  1.162      2.082
# mlr3::lrn('svm')                      1.172      2.063
# ---
#   The best learner for BLP is mlr3::lrn('svm') with lambda = 1.172.
# The best learner for GATES and CLAN is mlr3::lrn('ranger', num.trees = 100) with lambda.bar = 2.0818.


### 4.2 BLP ----
results_BLP <- get_BLP(genML)
results_BLP # print method
# BLP generic targets
# ---
#        Estimate CI lower CI upper p value
# beta.1   0.9750   0.8941    1.056  <2e-16 ***
# beta.2   1.0198   0.9442    1.095  <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# ---
# Confidence level of confidence interval [CI lower, CI upper]: 90 %

# the true ATE (= 0.979) is contained in the 90% confidence bounds of beta.1. Moreover, we reject the null of no significance of beta.2 at any reasonable level, which is expected since there is substantial treatment effect heterogeneity.

svg(file = "inst/doc/repo_plots/BLP.svg")
plot(results_BLP) # plot method
dev.off()


### 4.3 GATES ----
results_GATES <- get_GATES(genML)
results_GATES # print method
# GATES generic targets
# ---
#                 Estimate CI lower CI upper  p value
#  gamma.1          -0.4518  -0.6519   -0.256 8.16e-06 ***
#  gamma.2           0.4986   0.2974    0.697 1.07e-06 ***
#  gamma.3           1.4294   1.2280    1.631  < 2e-16 ***
#  gamma.4           2.4103   2.2066    2.612  < 2e-16 ***
#  gamma.4-gamma.1   2.8587   2.5754    3.142  < 2e-16 ***
#  ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  ---
#  Confidence level of confidence interval [CI lower, CI upper]: 90 %

plot(results_GATES) # plot method

svg(file = "inst/doc/repo_plots/GATES.svg")
plot(results_GATES) # plot method
dev.off()
# we know that here is treatment effect heterogeneity, so we expect a trend in the GATES per-group estimates as well as significance of all group differences. This is indeed the case.


### 4.4 CLAN of first variable ----
# plot HTEs against the variables
svg(file = "inst/doc/repo_plots/z1.svg")
plot(y = Z_CLAN[, "z1"], x = HTE, xlab = "True HTE", ylab = "Value of z1")
dev.off()

svg(file = "inst/doc/repo_plots/z2.svg")
plot(y = Z_CLAN[, "z2"], x = HTE, xlab = "True HTE", ylab = "Value of z2")
dev.off()

svg(file = "inst/doc/repo_plots/z3.svg")
plot(y = Z_CLAN[, "z3"], x = HTE, xlab = "True HTE", ylab = "Value of z3")
dev.off()

svg(file = "inst/doc/repo_plots/random.svg")
plot(y = Z_CLAN[, "random"], x = HTE, xlab = "True HTE", ylab = "Value of 'random'")
dev.off()

# CLAN
results_CLAN_z1 <- get_CLAN(genML, variable = "z1")
svg(file = "inst/doc/repo_plots/CLAN_z1.svg")
plot(results_CLAN_z1)
dev.off()

results_CLAN_z2 <- get_CLAN(genML, variable = "z2")
svg(file = "inst/doc/repo_plots/CLAN_z2.svg")
plot(results_CLAN_z2)
dev.off()

results_CLAN_z3 <- get_CLAN(genML, variable = "z3")
svg(file = "inst/doc/repo_plots/CLAN_z3.svg")
plot(results_CLAN_z3)
dev.off()

results_CLAN_random <- get_CLAN(genML, variable = "random")
svg(file = "inst/doc/repo_plots/CLAN_random.svg")
plot(results_CLAN_random)
dev.off()
