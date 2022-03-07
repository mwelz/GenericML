rm(list = ls()) ; cat("\014")

# load data
load(paste0(getwd(), "/inst/doc/repo_plots/plot_data.Rdata"))

### 4.1 analyze results ----
## print
x

## the line below returns the medians of the estimated  \Lambda and \bar{\Lambda}
x$best$overview

## get best learner for BLP
x$best$BLP
# "mlr3::lrn('svm')"

## get best learner for GATES and CLAN (this is the same learner)
x$best$GATES
x$best$CLAN
# "mlr3::lrn('ranger', num.trees = 100)"


### 4.2 BLP ----
get_BLP(x, plot = FALSE)
#          Estimate  CB lower CB upper      Pr(>|z|)
#  beta.1 0.9749574 0.8938969 1.056229 7.139626e-123
#  beta.2 1.0198128 0.9441887 1.095462 2.449502e-151

# the true ATE (= 0.979) is contained in the 90% confidence bounds of beta.1. Moreover, we reject the null of no significance of beta.2 at any reasonable level, which is expected since there is substantial treatment effect heterogeneity.

svg(file = paste0(getwd(), "/inst/doc/repo_plots/BLP.svg"))
plot(x, type = "BLP") # plot.GenericML() method
dev.off()


### 4.3 GATES ----
get_GATES(x, plot = FALSE)
svg(file = paste0(getwd(), "/inst/doc/repo_plots/GATES.svg"))
plot(x, type = "GATES")
dev.off()
# we know that here is treatment effect heterogeneity, so we expect a trend in the GATES per-group estimates as well as significance of all group differences. This is indeed the case.


### 4.4 CLAN of first variable ----
svg(file = paste0(getwd(), "/inst/doc/repo_plots/z1.svg"))
plot(y = Z_CLAN[, "z1"], x = HTE, xlab = "True HTE", ylab = "Value of z1")
dev.off()

svg(file = paste0(getwd(), "/inst/doc/repo_plots/z2.svg"))
plot(y = Z_CLAN[, "z2"], x = HTE, xlab = "True HTE", ylab = "Value of z2")
dev.off()

svg(file = paste0(getwd(), "/inst/doc/repo_plots/z3.svg"))
plot(y = Z_CLAN[, "z3"], x = HTE, xlab = "True HTE", ylab = "Value of z3")
dev.off()

svg(file = paste0(getwd(), "/inst/doc/repo_plots/random.svg"))
plot(y = Z_CLAN[, "random"], x = HTE, xlab = "True HTE", ylab = "Value of 'random'")
dev.off()

# CLAN
svg(file = paste0(getwd(), "/inst/doc/repo_plots/CLAN_z1.svg"))
get_CLAN(x, plot = TRUE, variable = "z1")
dev.off()

svg(file = paste0(getwd(), "/inst/doc/repo_plots/CLAN_z2.svg"))
get_CLAN(x, plot = TRUE, variable = "z2")
dev.off()

svg(file = paste0(getwd(), "/inst/doc/repo_plots/CLAN_z3.svg"))
get_CLAN(x, plot = TRUE, variable = "z3")
dev.off()

svg(file = paste0(getwd(), "/inst/doc/repo_plots/CLAN_random.svg"))
get_CLAN(x, plot = TRUE, variable = "random")
dev.off()
